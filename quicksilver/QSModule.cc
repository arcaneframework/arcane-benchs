// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "QSModule.hh"

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <sched.h>
#include "QS_Vector.hh"
#include "utilsMpi.hh"
#include "MonteCarlo.hh"
#include "MC_Processor_Info.hh"
#include "DecompositionObject.hh"
#include "GlobalFccGrid.hh"
#include "MeshPartition.hh"
#include "CommObject.hh"
#include "SharedMemoryCommObject.hh"
#include "MpiCommObject.hh"
#include "MC_Vector.hh"
#include "NuclearData.hh"
#include "MaterialDatabase.hh"
#include "MC_Time_Info.hh"
#include "Tallies.hh"
#include "MC_Base_Particle.hh"
#include "PhysicalConstants.hh"
#include "MCT.hh"





#include "arcane/IParallelMng.h"
#include "arcane/IMesh.h"

using namespace Arcane;
using namespace std;

#define MAX_PRODUCTION_SIZE 4

bool ordre_qs[]  = {true, true, false, false, false, true};

int QS2ArcaneFacet[] = {16, 17, 18, 19, 
                        7 , 6 , 5 , 4 , 
                        20, 23, 22, 21, 
                        8 , 9 , 10, 11, 
                        12, 13, 14, 15, 
                        3 , 2 , 1 , 0 };

int QS2ArcaneFace[] = {4, 1, 5, 2, 3, 0};

int QS2ArcaneNode[] = {0, 1, 2, 3,
                       0, 3, 2, 1,
                       1, 0, 3, 2,
                       0, 1, 2, 3,
                       0, 1, 2, 3,
                       0, 3, 2, 1};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void QSModule::
startInit()
{
  info() << "Module Quicksilver INIT"; 

  cartesian_mesh = ICartesianMesh::getReference(mesh(), true);
  // mpiInit(&argc, &argv);
  printBanner(GIT_VERS, GIT_HASH);
  int argc = 3;
  char *argv[] = {".", "-i", "/home/lheritiera/Documents/arcane/arcane-benchs/quicksilver/Coral2_P1.inp"};

  params = getParameters(argc, argv);
  cartesian_mesh->computeDirections();
  getParametersAxl();
  printParameters(params, cout);

  // monteCarlo stores just about everything.
  monteCarlo = initMC(params);
  monteCarloArc = initMCArc(params);

  MC_FASTTIMER_START(MC_Fast_Timer::main); // this can be done once monteCarlo exist.
}


void QSModule::
cycleInit()
{
  info() << "Module Quicksilver cycleInit";

  bool loadBalance = (bool)params.simulationParams.loadBalance;

  MC_FASTTIMER_START(MC_Fast_Timer::cycleInit);

  monteCarlo->clearCrossSectionCache();
  clearCrossSectionCache(); // Arc

  monteCarlo->_tallies->CycleInitialize(monteCarlo); // Ne fait rien.

  monteCarlo->_particleVaultContainer->swapProcessingProcessedVaults();
  monteCarlo->_particleVaultContainer->collapseProcessed();
  monteCarlo->_particleVaultContainer->collapseProcessing();

  // Nombre de particles dans Processing.
  monteCarlo->_tallies->_balanceTask[0]._start =
      monteCarlo->_particleVaultContainer->sizeProcessing();

  monteCarlo->particle_buffer->Initialize();

  MC_SourceNow(monteCarlo);
  MC_SourceNowArc(monteCarloArc);

  // Réduction ou augmentation du nombre de particules.
  PopulationControl(monteCarlo, loadBalance); // controls particle population
  PopulationControlArc(monteCarloArc, loadBalance); // controls particle population

  // Roulette sur les particules avec faible poids.
  RouletteLowWeightParticles(monteCarlo); // Delete particles with low statistical weight
  RouletteLowWeightParticlesArc(monteCarloArc); // Delete particles with low statistical weight

  MC_FASTTIMER_STOP(MC_Fast_Timer::cycleInit);
}

void QSModule::
cycleTracking()
{
  info() << "Module Quicksilver cycleTracking";

  MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking);

  //tracking(monteCarlo);

  info() << "Module Quicksilver cycleTrackingArc";

  trackingArc(monteCarloArc);

  MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking);
}

void QSModule::
cycleFinalize()
{
  info() << "Module Quicksilver cycleFinalize";

  MC_FASTTIMER_START(MC_Fast_Timer::cycleFinalize);

  monteCarlo->_tallies->_balanceTask[0]._end =
      monteCarlo->_particleVaultContainer->sizeProcessed();

  // Update the cumulative tally data.
  monteCarlo->_tallies->CycleFinalize(monteCarlo);

  monteCarlo->time_info->cycle++;

  monteCarlo->particle_buffer->Free_Memory();

  MC_FASTTIMER_STOP(MC_Fast_Timer::cycleFinalize);

  monteCarlo->fast_timer->Last_Cycle_Report(params.simulationParams.cycleTimers,
                                        monteCarlo->processor_info->rank,
                                        monteCarlo->processor_info->num_processors,
                                        monteCarlo->processor_info->comm_mc_world);
}

void QSModule::
gameOver()
{
  MC_FASTTIMER_STOP(MC_Fast_Timer::main);

  info() << "Module Quicksilver gameOver";

  monteCarlo->fast_timer->Cumulative_Report(
      monteCarlo->processor_info->rank, monteCarlo->processor_info->num_processors,
      monteCarlo->processor_info->comm_mc_world,
      monteCarlo->_tallies->_balanceCumulative._numSegments);
  monteCarlo->_tallies->_spectrum.PrintSpectrum(monteCarlo);

  coralBenchmarkCorrectness(monteCarlo, params);

  #ifdef HAVE_UVM
  monteCarlo->~MonteCarlo();
  cudaFree(monteCarlo);
#else
  delete monteCarlo;
#endif

  //mpiFinalize();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void QSModule::
getParametersAxl()
{
  // Equivalent de Parameters::parseCommandLine().
  params.simulationParams.dt = options()->getDt();
  params.simulationParams.fMax = options()->getFMax();
  params.simulationParams.nParticles = options()->getNParticles();
  params.simulationParams.nSteps = options()->getNSteps();
  params.simulationParams.seed = options()->getSeed();
  {
    CellDirectionMng cdm(cartesian_mesh->cellDirection(MD_DirX));
    params.simulationParams.nx = cdm.globalNbCell();
  }
  {
    CellDirectionMng cdm(cartesian_mesh->cellDirection(MD_DirY));
    params.simulationParams.ny = cdm.globalNbCell();
  }
  {
    CellDirectionMng cdm(cartesian_mesh->cellDirection(MD_DirZ));
    params.simulationParams.nz = cdm.globalNbCell();
  }

  // TODO : lx, ly, lz à calculer avec les valeurs de mesh.
  params.simulationParams.lx = options()->getLx();
  params.simulationParams.ly = options()->getLy();
  params.simulationParams.lz = options()->getLz();

  // xDom, yDom, zDom à calculer avec les valeurs de nb-part-x.

  // ???
  // addArg("bTally",           'B', 1, 'i', &(sp.balanceTallyReplications), 0, "number of balance tally replications");
  // addArg("fTally",           'F', 1, 'i', &(sp.fluxTallyReplications),    0, "number of scalar flux tally replications");
  // addArg("cTally",           'C', 1, 'i', &(sp.cellTallyReplications),    0, "number of scalar cell tally replications");



  // Equivalent de Parameters::scanSimulationBlock().
  params.simulationParams.boundaryCondition = options()->getBoundaryCondition().localstr();
  params.simulationParams.eMin = options()->getEMin();
  params.simulationParams.eMax = options()->getEMax();
  params.simulationParams.nGroups = options()->getNGroups();
  params.simulationParams.lowWeightCutoff = options()->getLowWeightCutoff();

}


MonteCarlo* QSModule::
initMCArc(const Parameters& params)
{
   MonteCarlo* monteCarlo;
   monteCarlo = new MonteCarlo(params);

   monteCarlo->time_info->time_step = params.simulationParams.dt;

   initNuclearData(monteCarlo, params); // Configuration des materiaux.
   initMesh(monteCarlo, params);
   initTallies(monteCarlo, params);

   MC_Base_Particle::Update_Counts();

   //   used when debugging cross sections
   checkCrossSections(monteCarlo, params);
   return monteCarlo;
}


/// Initializes both the NuclearData and the MaterialDatabase.  These
/// two structures are inherently linked since the isotopeGids stored in
/// the MaterialDatabase must correspond to the isotope indices in the
/// NuclearData.

void QSModule::
initNuclearData(MonteCarlo* monteCarlo, const Parameters& params)
{

  monteCarlo->_nuclearData = new NuclearData(params.simulationParams.nGroups,
                                                params.simulationParams.eMin,
                                                params.simulationParams.eMax);
  monteCarlo->_materialDatabase = new MaterialDatabase();

  map<string, Polynomial> crossSection;
  for (auto crossSectionIter = params.crossSectionParams.begin();
      crossSectionIter != params.crossSectionParams.end();
      crossSectionIter++)
  {
    const CrossSectionParameters& cp = crossSectionIter->second;
    crossSection.insert(make_pair(cp.name, Polynomial(cp.aa, cp.bb, cp.cc, cp.dd, cp.ee)));
  }
  
  int num_isotopes  = 0;
  int num_materials = 0;
  
  for( auto matIter = params.materialParams.begin(); matIter != params.materialParams.end(); matIter++ )
  {
    const MaterialParameters& mp = matIter->second;
    num_isotopes += mp.nIsotopes;
    num_materials++;
  }
  
  monteCarlo->_nuclearData->_isotopes.reserve( num_isotopes, VAR_MEM );
  monteCarlo->_materialDatabase->_mat.reserve( num_materials, VAR_MEM );
  
  for (auto matIter = params.materialParams.begin();
      matIter != params.materialParams.end(); matIter++)
  {
    const MaterialParameters& mp = matIter->second;
    Material material(mp.name, mp.mass);
    double nuBar = params.crossSectionParams.at(mp.fissionCrossSection).nuBar;
    material._iso.reserve( mp.nIsotopes, VAR_MEM );
    
    for (int iIso=0; iIso<mp.nIsotopes; ++iIso)
    {
        int isotopeGid = monteCarlo->_nuclearData->addIsotope(
          mp.nReactions,
          crossSection.at(mp.fissionCrossSection),
          crossSection.at(mp.scatteringCrossSection),
          crossSection.at(mp.absorptionCrossSection),
          nuBar,
          mp.totalCrossSection,
          mp.fissionCrossSectionRatio,
          mp.scatteringCrossSectionRatio,
          mp.absorptionCrossSectionRatio);
        
        // atomFraction for each isotope is 1/nIsotopes.  Treats all
        // isotopes as equally prevalent.
        material.addIsotope(Isotope(isotopeGid, 1.0/mp.nIsotopes));
    }
    monteCarlo->_materialDatabase->addMaterial(material);
  }
}




void QSModule::
initMesh(MonteCarlo* monteCarlo, const Parameters& params)
{
    
  int nx = params.simulationParams.nx;
  int ny = params.simulationParams.ny;
  int nz = params.simulationParams.nz;

  double lx = params.simulationParams.lx;
  double ly = params.simulationParams.ly;
  double lz = params.simulationParams.lz;

  int xDom = params.simulationParams.xDom;
  int yDom = params.simulationParams.yDom;
  int zDom = params.simulationParams.zDom;

  double dx = lx / nx;
  double dy = ly / ny;
  double dz = lz / nz;

  int myRank, nRanks;
  myRank = meshHandle().subDomain()->parallelMng()->commRank();
  nRanks = meshHandle().subDomain()->parallelMng()->commSize();

  // Extrait de GlobalFccGrid.cc.
  static vector<Tuple4> offset;
  offset.reserve(14);
  offset.push_back(Tuple4(0, 0, 0, 0)); // 0
  offset.push_back(Tuple4(1, 0, 0, 0)); // 1
  offset.push_back(Tuple4(1, 1, 0, 0)); // 3
  offset.push_back(Tuple4(0, 1, 0, 0)); // 2

  offset.push_back(Tuple4(0, 0, 1, 0)); // 4
  offset.push_back(Tuple4(1, 0, 1, 0)); // 5
  offset.push_back(Tuple4(1, 1, 1, 0)); // 7
  offset.push_back(Tuple4(0, 1, 1, 0)); // 6

  offset.push_back(Tuple4(0, 0, 0, 3)); // 13
  offset.push_back(Tuple4(0, 0, 0, 1)); // 9
  offset.push_back(Tuple4(0, 0, 0, 2)); // 11
  offset.push_back(Tuple4(0, 0, 1, 3)); // 12
  offset.push_back(Tuple4(1, 0, 0, 1)); // 8
  offset.push_back(Tuple4(0, 1, 0, 2)); // 10

  static vector<Tuple4> faceTupleOffset;
  faceTupleOffset.reserve(6);
  faceTupleOffset.push_back( Tuple4( 0,  0, -1, 5) ); // 13
  faceTupleOffset.push_back( Tuple4(-1,  0,  0, 1) ); // 9
  faceTupleOffset.push_back( Tuple4( 0, -1,  0, 3) ); // 11
  faceTupleOffset.push_back( Tuple4( 0,  0,  1, 4) ); // 12
  faceTupleOffset.push_back( Tuple4( 1,  0,  0, 0) ); // 8
  faceTupleOffset.push_back( Tuple4( 0,  1,  0, 2) ); // 10
   

  // Ici, les cells ont déjà une numérotation (cell.uniqueId() (G)  ou  cell.localId() (L))
  // On doit numéroter les nodes selon cell.uniqueId.
  //////////////////// Début Numérotation Node/Face /////////////////////////
  #define VERIF true
  #if VERIF
    m_coord.resize(5);
  #else
    m_coord.resize(4);
  #endif
  m_coordCm.resize(3);
  m_coordMidCm.resize(3);
  m_coordMid.resize(4);
  m_coordFace.resize(4);
  m_coordCenter.resize(4);
  m_total.resize(monteCarlo->_nuclearData->_numEnergyGroups);
  ENUMERATE_CELL(icell, ownCells())
  {
    Cell cell = *icell;
    Int32 unique_id = cell.uniqueId().asInt32();

    int index = unique_id;
    int x = index % nx;
    index /= nx;
    int y = index % ny;
    int z = index / ny;
    int compt_face = 13;

    m_coordCenter[icell][0] = 0.0;
    m_coordCenter[icell][1] = 0.0;
    m_coordCenter[icell][2] = 0.0;

    int compt = 0;
    ENUMERATE_NODE(inode, cell.nodes())
    {
      #if VERIF
      if(m_coord[inode][4] == -123)
      {
        if( m_coord[inode][0] != offset[compt].x() + x ||
            m_coord[inode][1] != offset[compt].y() + y ||
            m_coord[inode][2] != offset[compt].z() + z ||
            m_coord[inode][3] != offset[compt].b())
          {
            error() << m_coord[inode][0] << " " << offset[compt].x() + x;
            error() << m_coord[inode][1] << " " << offset[compt].y() + y;
            error() << m_coord[inode][2] << " " << offset[compt].z() + z;
            error() << m_coord[inode][3] << " " << offset[compt].b();
            ARCANE_FATAL("Erreur Node");
          }
      }
      m_coord[inode][4] = -123;
      #endif
      m_coord[inode][0] = offset[compt].x() + x;
      m_coord[inode][1] = offset[compt].y() + y;
      m_coord[inode][2] = offset[compt].z() + z;
      m_coord[inode][3] = offset[compt].b();
      compt++;
      //info() << m_coord[inode][0] << " " << m_coord[inode][1] << " " << m_coord[inode][2];

      m_coordCm[inode][0] = m_coord[inode][0] * dx;
      m_coordCm[inode][1] = m_coord[inode][1] * dy;
      m_coordCm[inode][2] = m_coord[inode][2] * dz;

      m_coordCenter[icell][0] += m_coordCm[inode][0];
      m_coordCenter[icell][1] += m_coordCm[inode][1];
      m_coordCenter[icell][2] += m_coordCm[inode][2];
    }
    m_coordCenter[icell][0] /= cell.nbNode();
    m_coordCenter[icell][1] /= cell.nbNode();
    m_coordCenter[icell][2] /= cell.nbNode();

    double volume = 0;
    MC_Vector cellCenter(m_coordCenter[icell][0], m_coordCenter[icell][1], m_coordCenter[icell][2]);

    ENUMERATE_FACE(iface, cell.faces())
    {
      Face face = *iface;

      m_indexArc[iface] = iface.index();

      m_coordMid[iface][0] = offset[compt].x() + x;
      m_coordMid[iface][1] = offset[compt].y() + y;
      m_coordMid[iface][2] = offset[compt].z() + z;
      m_coordMid[iface][3] = offset[compt].b();

      m_coordMidCm[iface][0] = m_coordMid[iface][0] * dx;
      m_coordMidCm[iface][1] = m_coordMid[iface][1] * dy;
      m_coordMidCm[iface][2] = m_coordMid[iface][2] * dz;

      if(m_coordMid[iface][3] == 1)
      {
        m_coordMidCm[iface][1] += dy / 2;
        m_coordMidCm[iface][2] += dz / 2;
      }
      else if(m_coordMid[iface][3] == 2)
      {
        m_coordMidCm[iface][0] += dx / 2;
        m_coordMidCm[iface][2] += dz / 2;
      }
      else
      {
        m_coordMidCm[iface][0] += dx / 2;
        m_coordMidCm[iface][1] += dy / 2;
      }

      int compt2 = compt - 8;

      m_coordFace[iface][0] = std::min(std::max(0, faceTupleOffset[compt2].x() + x), nx-1); // MinMax pour éviter pos négatives.
      m_coordFace[iface][1] = std::min(std::max(0, faceTupleOffset[compt2].y() + y), ny-1);
      m_coordFace[iface][2] = std::min(std::max(0, faceTupleOffset[compt2].z() + z), nz-1);


      // Définir les conditions boundary.
      //////////////////// Début conditions boundary /////////////////////////
      qs_vector<MC_Subfacet_Adjacency_Event::Enum> condiBound = getBoundaryCondition(params);
      //Real faceNbr = m_coordFace[iface][0] + nx*(m_coordFace[iface][1] + ny*(m_coordFace[iface][2]));

      // Si la face est au bord du domaine entier.
      if(m_coordFace[iface][0] == x && m_coordFace[iface][1] == y && m_coordFace[iface][2] == z)
      {
        m_boundaryCond[iface] = condiBound[faceTupleOffset[compt2].b()];
      }
      else
      {
        Face face = *iface;
        // Si la face est au bord du sous-domaine.
        if(face.isSubDomainBoundary())
        {
          info() <<  "x : " << x << " y: " << y << " z: " << z << " xx : " << faceTupleOffset[compt2].x() << " yy : "<< faceTupleOffset[compt2].y() << " zz : " << faceTupleOffset[compt2].z() << " compt : "<< compt2;
          ARCANE_FATAL("TODO Impossible");
          m_boundaryCond[iface] = MC_Subfacet_Adjacency_Event::Transit_Off_Processor;
        }
        
        else
        {
          m_boundaryCond[iface] = MC_Subfacet_Adjacency_Event::Transit_On_Processor;
        }
      }
      //////////////////// Fin conditions boundary /////////////////////////

      //////////////////// Début Volume Cell /////////////////////////


      for (int i = 0; i < 4; i++)
      {
        int first_pos_node = (ordre_qs[iface.index()] ? ((i == 3) ? 0 : i+1) : i);
        int second_pos_node = (ordre_qs[iface.index()] ? i : ((i == 3) ? 0 : i+1));

        Node first_node = face.node(first_pos_node);
        Node second_node = face.node(second_pos_node);

        MC_Vector aa = MC_Vector(m_coordCm[first_node][0], m_coordCm[first_node][1], m_coordCm[first_node][2]) - cellCenter;
        MC_Vector bb = MC_Vector(m_coordCm[second_node][0], m_coordCm[second_node][1], m_coordCm[second_node][2]) - cellCenter;
        MC_Vector cc = MC_Vector(m_coordMidCm[iface][0], m_coordMidCm[iface][1], m_coordMidCm[iface][2]) - cellCenter;

        volume += abs(aa.Dot(bb.Cross(cc)));
      }

      
      compt++;
    }
    volume /= 6.0;

    m_volume[cell] = volume;
    //////////////////// Fin Volume Cell /////////////////////////


    std::string matName = findMaterial(params, cellCenter);
    m_material[icell] = monteCarlo->_materialDatabase->findMaterial(matName);

    for(int i = 0; i < monteCarlo->_nuclearData->_numEnergyGroups; i++)
    {
      m_total[icell][i] = 0.0;
    }

    m_cellNumberDensity[icell] = 1.0;
    m_sourceTally[icell] = 0;

    
  }
  //////////////////// Fin Numérotation Node/Face /////////////////////////
  
  if (myRank == 0) { cout << "Finished initMesh" <<endl; }
}

qs_vector<MC_Subfacet_Adjacency_Event::Enum> QSModule::
getBoundaryCondition(const Parameters& params)
   {
      qs_vector<MC_Subfacet_Adjacency_Event::Enum> bc(6);
      if (params.simulationParams.boundaryCondition == "reflect")
         bc = qs_vector<MC_Subfacet_Adjacency_Event::Enum>(6, MC_Subfacet_Adjacency_Event::Boundary_Reflection);
      else if (params.simulationParams.boundaryCondition == "escape")
         bc = qs_vector<MC_Subfacet_Adjacency_Event::Enum>(6, MC_Subfacet_Adjacency_Event::Boundary_Escape);
      else if  (params.simulationParams.boundaryCondition == "octant")
         for (unsigned ii=0; ii<6; ++ii)
         {
            if (ii % 2 == 0) bc[ii] = MC_Subfacet_Adjacency_Event::Boundary_Escape;
            if (ii % 2 == 1) bc[ii] = MC_Subfacet_Adjacency_Event::Boundary_Reflection;
         }
      else
         qs_assert(false);
      return bc;
   }



void QSModule::
initTallies(MonteCarlo* monteCarlo, const Parameters& params)
{
  // Pour l'instant, balanceTallyReplications = fluxTallyReplications = cellTallyReplications = 1.
  m_absorb = 0.0;      // Number of particles absorbed
  m_census = 0.0;      // Number of particles that enter census
  m_escape = 0.0;      // Number of particles that escape
  m_collision = 0.0;   // Number of collosions
  m_end = 0.0;         // Number of particles at end of cycle
  m_fission = 0.0;     // Number of fission events
  m_produce = 0.0;     // Number of particles created by collisions
  m_scatter = 0.0;     // Number of scatters
  m_start = 0.0;       // Number of particles at beginning of cycle
  m_source = 0.0;      // Number of particles sourced in
  m_rr = 0.0;          // Number of particles Russian Rouletted in population control
  m_split = 0.0;       // Number of particles split in population control
  m_numSegments = 0.0; // Number of segements

  int sizeOfSFT = monteCarlo->_nuclearData->_energies.size()-1;
  m_scalarFluxTally.resize(sizeOfSFT);

  ENUMERATE_CELL(icell, ownCells())
  {
    m_cellTally[icell] = 0.0;
    for (int i = 0; i < sizeOfSFT; i++)
    {
      m_scalarFluxTally[icell][i] = 0.0;
    }
  }


}




void QSModule::
consistencyCheck(int myRank, const qs_vector<MC_Domain>& domain)
{
  if (myRank == 0) { cout << "Starting Consistency Check" <<endl; }
  unsigned nDomains = domain.size();
  for (int iDomain=0; iDomain<nDomains; ++iDomain)
  {
      const MC_Mesh_Domain& mesh = domain[iDomain].mesh;
      unsigned nCells = mesh._cellConnectivity.size();
      for (unsigned iCell=0; iCell<nCells; ++iCell)
      {
        for (unsigned iFacet = 0; iFacet<24; ++iFacet)
        {
            const MC_Location& current =
              mesh._cellConnectivity[iCell]._facet[iFacet].subfacet.current;
            qs_assert(current.cell == iCell);

            const MC_Location& adjacent =
              mesh._cellConnectivity[iCell]._facet[iFacet].subfacet.adjacent;

            int jDomain = adjacent.domain;
            int jCell = adjacent.cell;
            int jFacet = adjacent.facet;

            const Subfacet_Adjacency& backside = domain[jDomain].mesh._cellConnectivity[jCell]._facet[jFacet].subfacet;

            qs_assert (backside.adjacent.domain == iDomain);
            qs_assert (backside.adjacent.cell == iCell);
            qs_assert (backside.adjacent.facet == iFacet);
        }
      }
  }
  if (myRank == 0) { cout << "Finished Consistency Check" <<endl; }
}




// scatter the centers (somewhat) randomly
void QSModule::
initializeCentersRandomly(int nCenters,
                      const GlobalFccGrid& grid,
                      vector<MC_Vector>& centers)
{
  set<Tuple> picked;
  do
  {
      Tuple iTuple(drand48()*grid.nx()/2,
                  drand48()*grid.ny()/2,
                  drand48()*grid.nz()/2);

      if (!picked.insert(iTuple).second)
        continue;

      iTuple += iTuple; // iTuple *= 2;
      Long64 iCell = grid.cellTupleToIndex(iTuple);
      MC_Vector r = grid.cellCenter(iCell);
      centers.push_back(r);
  } while (centers.size() < nCenters);
}



void QSModule::
initializeCentersGrid(double lx, double ly, double lz,
                          int xDom, int yDom, int zDom,
                          vector<MC_Vector>& centers)
{
  double dx = lx/xDom;
  double dy = ly/yDom;
  double dz = lz/zDom;
  for (int ix=0; ix<xDom; ++ix)
      for (int iy=0; iy<yDom; ++iy)
        for (int iz=0; iz<zDom; ++iz)
            centers.push_back(
              MC_Vector( (0.5+ix)*dx, (0.5+iy)*dy, (0.5+iz)*dz )
            );
}

// Returns true if the specified coordinate in inside the specified
// geometry.  False otherwise
bool QSModule::
isInside(const GeometryParameters& geom, const MC_Vector& rr)
{
  bool inside = false;
  switch (geom.shape)
  {
    case GeometryParameters::BRICK:
      {
        if ( (rr.x >= geom.xMin && rr.x <= geom.xMax) &&
              (rr.y >= geom.yMin && rr.y <= geom.yMax) &&
              (rr.z >= geom.zMin && rr.z <= geom.zMax) )
            inside = true;
      }
      break;
    case GeometryParameters::SPHERE:
      {
        MC_Vector center(geom.xCenter, geom.yCenter, geom.zCenter);
        if ( (rr-center).Length() <= geom.radius)
            inside = true;
      }

      break;
    default:
      qs_assert(false);
  }
  return inside;
}

// Returns the name of the material present at coordinate rr.  If
// multiple materials overlap return the last material found.
string QSModule::
findMaterial(const Parameters& params, const MC_Vector& rr)
{
  string materialName;
  for (unsigned ii=0; ii< params.geometryParams.size(); ++ii)
      if (isInside(params.geometryParams[ii], rr))
        materialName = params.geometryParams[ii].materialName;

  qs_assert(materialName.size() > 0);
  return materialName;
}

// This function is useful for debugging but is not called in ordinary
// use of the code.  Uncomment the call to this function in initMC()
// if you want to get plot data for the cross sections.
void QSModule::
checkCrossSections(MonteCarlo* monteCarlo, const Parameters& params)
{
  if( monteCarlo->_params.simulationParams.crossSectionsOut == "" ) return;

  struct XC_Data
  {
      XC_Data() : absorption(0.), fission(0.), scatter(0.){}
      double absorption;
      double fission;
      double scatter;
  };

  NuclearData* nd = monteCarlo->_nuclearData;
  int nGroups = nd->_energies.size() - 1;
  vector<double> energy(nGroups);
  for (unsigned ii=0; ii<nGroups; ++ii)
      energy[ii] = (nd->_energies[ii] + nd->_energies[ii+1])/2.0;


  MaterialDatabase* matDB = monteCarlo->_materialDatabase;
  unsigned nMaterials = matDB->_mat.size();

  map<string, vector<XC_Data> > xcTable;


  // for each material
  for (unsigned iMat=0; iMat<nMaterials; ++iMat)
  {
      const string& materialName = matDB->_mat[iMat]._name;
      vector<XC_Data>& xcVec = xcTable[materialName];
      xcVec.resize(nGroups);
      unsigned nIsotopes = matDB->_mat[iMat]._iso.size();
      // for each isotope
      for (unsigned iIso=0; iIso<nIsotopes; ++iIso)
      {
        int isotopeGid = monteCarlo->_materialDatabase->_mat[iMat]._iso[iIso]._gid;
        unsigned nReactions = nd->_isotopes[isotopeGid]._species[0]._reactions.size();
        // for each reaction
        for (unsigned iReact=0; iReact<nReactions; ++iReact)
        {
            // loop over energies
            NuclearDataReaction& reaction = nd->_isotopes[isotopeGid]._species[0]._reactions[iReact];
            // accumulate cross sections by reaction type
            for (unsigned iGroup=0; iGroup<nGroups; ++iGroup)
            {
              switch (reaction._reactionType)
              {
                case NuclearDataReaction::Scatter:
                  xcVec[iGroup].scatter += reaction.getCrossSection(iGroup)/nIsotopes;
                  break;
                case NuclearDataReaction::Absorption:
                  xcVec[iGroup].absorption += reaction.getCrossSection(iGroup)/nIsotopes;
                  break;
                case NuclearDataReaction::Fission:
                  xcVec[iGroup].fission += reaction.getCrossSection(iGroup)/nIsotopes;
                  break;
                case NuclearDataReaction::Undefined:
                  qs_assert(false);
                  break;
              }   
            }
        }
      }
  }

FILE* xSec;

std::string fileName = monteCarlo->_params.simulationParams.crossSectionsOut + ".dat";

xSec = fopen( fileName.c_str(), "w" );

  // print cross section data
  // first the header
  fprintf(xSec, "#group  energy");
  for (auto mapIter=xcTable.begin(); mapIter!=xcTable.end(); ++mapIter)
  {
      const string& materialName = mapIter->first;
      fprintf(xSec, "  %s_a  %s_f  %s_s", materialName.c_str(), materialName.c_str(), materialName.c_str());
  }
  fprintf(xSec,"\n");

  // now the data
  for (unsigned ii=0; ii<nGroups; ++ii)
  {
      fprintf(xSec, "%u  %g", ii, energy[ii]);
      for (auto mapIter=xcTable.begin(); mapIter!=xcTable.end(); ++mapIter)
      {
        fprintf(xSec, "  %g  %g  %g", mapIter->second[ii].absorption, mapIter->second[ii].fission, mapIter->second[ii].scatter);
      }
      fprintf(xSec, "\n");
  }
fclose( xSec );
}

void QSModule::
clearCrossSectionCache()
{
  ENUMERATE_CELL(icell, ownCells())
  {
    for(int i = 0; i < monteCarlo->_nuclearData->_numEnergyGroups; i++)
    {
      m_total[icell][i] = 0.0;
    }
  }
}

void QSModule::
MC_SourceNowArc(MonteCarlo *monteCarlo)
{
    NVTX_Range range("MC_Source_Now");
  
    std::vector<double> source_rate(monteCarlo->_materialDatabase->_mat.size());  // Get this from user input

    for ( int material_index = 0; material_index < monteCarlo->_materialDatabase->_mat.size(); material_index++ )
    {
        std::string name = monteCarlo->_materialDatabase->_mat[material_index]._name;
        double sourceRate = monteCarlo->_params.materialParams[name].sourceRate;
        source_rate[material_index] = sourceRate;
    }

    double local_weight_particles = 0;

    ENUMERATE_CELL(icell, ownCells())
    {
      double cell_weight_particles = m_volume[icell] * source_rate[m_material[icell]] * monteCarlo->time_info->time_step;
      local_weight_particles += cell_weight_particles;
    }

    double total_weight_particles = 0;

    total_weight_particles = mesh()->parallelMng()->reduce(Parallel::ReduceSum, local_weight_particles);

    uint64_t num_particles = monteCarlo->_params.simulationParams.nParticles;
    double source_fraction = 0.1;
    double source_particle_weight = total_weight_particles/(source_fraction * num_particles);
    // Store the source particle weight for later use.
    monteCarlo->source_particle_weight = source_particle_weight;

    uint64_t vault_size       = monteCarlo->_particleVaultContainer->getVaultSize();
    uint64_t processing_index = monteCarlo->_particleVaultContainer->sizeProcessing() / vault_size;

    uint64_t task_index = 0;
    uint64_t particle_count = 0;

    ENUMERATE_CELL(icell, ownCells())
    {
      Cell cell = *icell;
      double cell_weight_particles = m_volume[icell] * source_rate[m_material[icell]] * monteCarlo->time_info->time_step;
      double cell_num_particles_float = cell_weight_particles / source_particle_weight;
      int cell_num_particles = (int)cell_num_particles_float;

      for ( int particle_index = 0; particle_index < cell_num_particles; particle_index++ )
      {
        MC_Particle particle;
        uint64_t random_number_seed;

        // TODO : Voir si c'est utile de faire des atomics (non ?).
        ATOMIC_CAPTURE( m_sourceTally[icell], 1, random_number_seed );
        random_number_seed += cell.uniqueId().asInt64() * UINT64_C(0x0100000000); // Voir à quoi sert la mult (MC_Domain::MC_Domain).

        particle.random_number_seed = rngSpawn_Random_Number_Seed(&random_number_seed);
        particle.identifier = random_number_seed;

        // cout << "particle.identifier : " << particle.identifier << endl;

        MCT_Generate_Coordinate_3D_GArc(&particle.random_number_seed, cell, particle.coordinate, monteCarlo);

        // if(particle_index < 11)
        // {
        //   cout << "particle.identifier : " << particle.identifier << endl;
        //   cout << particle.coordinate.x << " x " << particle.coordinate.y << " x " << particle.coordinate.z << endl;
        // }
        // else
        // {
        //   exit(123);
        // }

        particle.direction_cosine.Sample_Isotropic(&particle.random_number_seed);
        particle.kinetic_energy = (monteCarlo->_params.simulationParams.eMax - monteCarlo->_params.simulationParams.eMin)*
                                rngSample(&particle.random_number_seed) + monteCarlo->_params.simulationParams.eMin;

        double speed = Get_Speed_From_Energy(particle.kinetic_energy);

        particle.velocity.x = speed * particle.direction_cosine.alpha;
        particle.velocity.y = speed * particle.direction_cosine.beta;
        particle.velocity.z = speed * particle.direction_cosine.gamma;

        particle.domain = mesh()->subDomain()->subDomainId();
        particle.cell   = cell.localId(); // TODO : Localid ou uniqueId ?
        particle.cellArc= cell;
        particle.task   = task_index;
        particle.weight = source_particle_weight;

        double randomNumber = rngSample(&particle.random_number_seed);
        particle.num_mean_free_paths = -1.0*std::log(randomNumber);

        randomNumber = rngSample(&particle.random_number_seed);
        particle.time_to_census = monteCarlo->time_info->time_step * randomNumber;

        MC_Base_Particle base_particle( particle );

        monteCarlo->_particleVaultContainer->addProcessingParticle( base_particle, processing_index );

        particle_count++;

        m_source = m_source() + 1; // TODO : Voir si besoin atomic.
      }
    }
}

double QSModule::
Get_Speed_From_Energy(double energy)
{
  static const double rest_mass_energy = PhysicalConstants::_neutronRestMassEnergy;
  static const double speed_of_light  = PhysicalConstants::_speedOfLight;


  return speed_of_light * sqrt(energy * (energy + 2.0*(rest_mass_energy)) /
                                ((energy + rest_mass_energy) * (energy + rest_mass_energy)));
}

void QSModule::
MCT_Generate_Coordinate_3D_GArc(uint64_t *random_number_seed,
                                  Cell &cell,
                                  MC_Vector &coordinate,  
                                  MonteCarlo* monteCarlo )
{
  // Determine the cell-center nodal point coordinates.
  MC_Vector center(m_coordCenter[cell][0], m_coordCenter[cell][1], m_coordCenter[cell][2]);

  int num_facets = 24;

  double random_number = rngSample(random_number_seed);
  double which_volume = random_number * 6.0 * m_volume[cell];

  // Find the tet to sample from.
  double current_volume = 0.0;
  int facet_index = -1;
  Node first_node;
  Node second_node;
  Face face;

  for(int i = 0; i < 6; i++)
  {
    face = cell.face(QS2ArcaneFace[i]);

    for (int j = 0; j < 4; j++)
    {
      facet_index++;

      int first_pos_node = QS2ArcaneNode[i*4 + j];
      int second_pos_node = QS2ArcaneNode[i*4 + ((j == 3) ? 0 : j+1)];

      first_node = face.node(first_pos_node);
      second_node = face.node(second_pos_node);

      MC_Vector point0(m_coordCm[first_node][0], m_coordCm[first_node][1], m_coordCm[first_node][2]);
      MC_Vector point1(m_coordCm[second_node][0], m_coordCm[second_node][1], m_coordCm[second_node][2]);
      MC_Vector point2(m_coordMidCm[face][0], m_coordMidCm[face][1], m_coordMidCm[face][2]);

      // cout << "i : " << i << " j : " << j << endl;
      // cout << "first_node[] : " << i*4 + j << " second_node[] : " << i*4 + ((j == 3) ? 0 : j+1) << endl;
      // cout << "first_node : " << first_pos_node << " second_node : " << second_pos_node << endl;
      // cout << point0.x << " x " << point0.y << " x " << point0.z << endl;
      // cout << point1.x << " x " << point1.y << " x " << point1.z << endl;
      // cout << point2.x << " x " << point2.y << " x " << point2.z << endl << endl;

      double subvolume = MCT_Cell_Volume_3D_G_vector_tetDetArc(point0, point1, point2, center);
      current_volume += subvolume;

      if(current_volume >= which_volume) { break; }
    }
    if(current_volume >= which_volume) { break; }
  }
  //cout << "current_volume : " << current_volume << endl;
  //exit(123);

  // Sample from the tet.
  double r1 = rngSample(random_number_seed);
  double r2 = rngSample(random_number_seed);
  double r3 = rngSample(random_number_seed);

  // Cut and fold cube into prism.
  if (r1 + r2 > 1.0)
  {
      r1 = 1.0 - r1;
      r2 = 1.0 - r2;
  }
  // Cut and fold prism into tetrahedron.
  if (r2 + r3 > 1.0)
  {
      double tmp = r3;
      r3 = 1.0 - r1 - r2;
      r2 = 1.0 - tmp;
  }
  else if (r1 + r2 + r3 > 1.0)
  {
      double tmp = r3;
      r3 = r1 + r2 + r3 - 1.0;
      r1 = 1.0 - r2 - tmp;
  }

  // numbers 1-4 are the barycentric coordinates of the random point.
  double r4 = 1.0 - r1 - r2 - r3;

  MC_Vector point0(m_coordCm[first_node][0], m_coordCm[first_node][1], m_coordCm[first_node][2]);
  MC_Vector point1(m_coordCm[second_node][0], m_coordCm[second_node][1], m_coordCm[second_node][2]);
  MC_Vector point2(m_coordMidCm[face][0], m_coordMidCm[face][1], m_coordMidCm[face][2]);

  // TODO : Pourrai ne pas être identique à QS_ori (pas le même ordre d'exploration).
  coordinate.x = ( r4 * center.x + r1 * point0.x + r2 * point1.x + r3 * point2.x );
  coordinate.y = ( r4 * center.y + r1 * point0.y + r2 * point1.y + r3 * point2.y );
  coordinate.z = ( r4 * center.z + r1 * point0.z + r2 * point1.z + r3 * point2.z );
}

///  \return 6 times the volume of the tet.
///
///  subtract v3 from v0, v1 and v2.  Then take the triple product of v0, v1 and v2.
double QSModule::
MCT_Cell_Volume_3D_G_vector_tetDetArc(const MC_Vector &v0_,
                                            const MC_Vector &v1_,
                                            const MC_Vector &v2_,
                                            const MC_Vector &v3)
{
  MC_Vector v0(v0_), v1(v1_), v2(v2_);

  v0.x -= v3.x; v0.y -= v3.y; v0.z -= v3.z;
  v1.x -= v3.x; v1.y -= v3.y; v1.z -= v3.z;
  v2.x -= v3.x; v2.y -= v3.y; v2.z -= v3.z;

  return
    v0.z*(v1.x*v2.y - v1.y*v2.x) +
    v0.y*(v1.z*v2.x - v1.x*v2.z) +
    v0.x*(v1.y*v2.z - v1.z*v2.y);
}

void QSModule::
PopulationControlArc(MonteCarlo* monteCarlo, bool loadBalance)
{
    NVTX_Range range("PopulationControl");

    uint64_t targetNumParticles = monteCarlo->_params.simulationParams.nParticles;
    uint64_t globalNumParticles = 0;
    uint64_t localNumParticles = monteCarlo->_particleVaultContainer->sizeProcessing();
   
    if (loadBalance)
    {
      // If we are parallel, we will have one domain per mpi processs.  The targetNumParticles is across
      // all MPI processes, so we need to divide by the number or ranks to get the per-mpi-process number targetNumParticles
      targetNumParticles = ceil((double)targetNumParticles / mesh()->parallelMng()->commSize() );

      //NO LONGER SPLITING VAULTS BY THREADS
//        // If we are threaded, targetNumParticles should be divided by the number of threads (tasks) to balance
//        // the particles across the thread level vaults.
//        targetNumParticles = ceil((double)targetNumParticles / (double)monteCarlo->processor_info->num_tasks);
    }
    else
    {
      globalNumParticles = mesh()->parallelMng()->reduce(Parallel::ReduceSum, localNumParticles);
    }
     
    Balance & taskBalance = monteCarlo->_tallies->_balanceTask[0];

    double splitRRFactor = 1.0;
    if (loadBalance)
    {
        int currentNumParticles = localNumParticles;
        if (currentNumParticles != 0)
            splitRRFactor = (double)targetNumParticles / (double)currentNumParticles;
        else
            splitRRFactor = 1.0;
    }
    else
    {
        splitRRFactor = (double)targetNumParticles / (double)globalNumParticles;
    }

    // On augmente ou diminue la population selon splitRRFactor (si > 1, on augmente en splittant ; si < 1, on diminue en killant ou en augmentant le poids (rand))
    if (splitRRFactor != 1.0)  // no need to split if population is already correct.
        PopulationControlGutsArc(splitRRFactor, localNumParticles, monteCarlo->_particleVaultContainer);

    monteCarlo->_particleVaultContainer->collapseProcessing();

    return;
}

void QSModule::
PopulationControlGutsArc(const double splitRRFactor, uint64_t currentNumParticles, ParticleVaultContainer* my_particle_vault)
{
    uint64_t vault_size = my_particle_vault->getVaultSize();
    uint64_t fill_vault_index = currentNumParticles / vault_size;

    // March backwards through the vault so killed particles doesn't mess up the indexing
    for (int particleIndex = currentNumParticles-1; particleIndex >= 0; particleIndex--)
    {
        uint64_t vault_index = particleIndex / vault_size; 

        ParticleVault& taskProcessingVault = *( my_particle_vault->getTaskProcessingVault(vault_index) );

        uint64_t taskParticleIndex = particleIndex%vault_size;

        MC_Base_Particle &currentParticle = taskProcessingVault[taskParticleIndex];
        double randomNumber = rngSample(&currentParticle.random_number_seed);
        if (splitRRFactor < 1)
        {
            if (randomNumber > splitRRFactor)
            {
                // Kill
	            taskProcessingVault.eraseSwapParticle(taskParticleIndex);
	            m_rr = m_rr() + 1; 
	        }
	        else
	        {
	            currentParticle.weight /= splitRRFactor;
	        }
        }
        else if (splitRRFactor > 1)
        {
            // Split
	        int splitFactor = (int)floor(splitRRFactor);
	        if (randomNumber > (splitRRFactor - splitFactor)) { splitFactor--; }
	  
	        currentParticle.weight /= splitRRFactor;
	        MC_Base_Particle splitParticle = currentParticle;
	  
	        for (int splitFactorIndex = 0; splitFactorIndex < splitFactor; splitFactorIndex++)
	        {
	            m_split = m_split() + 1;
	     
	            splitParticle.random_number_seed = rngSpawn_Random_Number_Seed(
			        &currentParticle.random_number_seed);
	            splitParticle.identifier = splitParticle.random_number_seed;

                my_particle_vault->addProcessingParticle( splitParticle, fill_vault_index );

	        }
        }
    }
}

void QSModule::
RouletteLowWeightParticlesArc(MonteCarlo* monteCarlo)
{
    NVTX_Range range("RouletteLowWeightParticles");

    const double lowWeightCutoff = monteCarlo->_params.simulationParams.lowWeightCutoff;

    if (lowWeightCutoff > 0.0)
    {

        uint64_t currentNumParticles = monteCarlo->_particleVaultContainer->sizeProcessing();
        uint64_t vault_size          = monteCarlo->_particleVaultContainer->getVaultSize();


	    // March backwards through the vault so killed particles don't mess up the indexing
	    const double source_particle_weight = monteCarlo->source_particle_weight;
	    const double weightCutoff = lowWeightCutoff*source_particle_weight;

	    for ( int64_t particleIndex = currentNumParticles-1; particleIndex >= 0; particleIndex--)
	    {
            uint64_t vault_index = particleIndex / vault_size; 

            ParticleVault& taskProcessingVault = *(monteCarlo->_particleVaultContainer->getTaskProcessingVault(vault_index));
            uint64_t taskParticleIndex = particleIndex%vault_size;
	        MC_Base_Particle &currentParticle = taskProcessingVault[taskParticleIndex];

	        if (currentParticle.weight <= weightCutoff)
	        {
	            double randomNumber = rngSample(&currentParticle.random_number_seed);
	            if (randomNumber <= lowWeightCutoff)
	            {
		            // The particle history continues with an increased weight.
		            currentParticle.weight /= lowWeightCutoff;
	            }
	            else
	            {
		            // Kill
		            taskProcessingVault.eraseSwapParticle(taskParticleIndex);
		            m_rr = m_rr() + 1;
	            } 
	        }
	    }
        monteCarlo->_particleVaultContainer->collapseProcessing();
    }
}

void QSModule::
tracking(MonteCarlo* monteCarlo)
{
  bool done = false;

  // Determine whether or not to use GPUs if they are available (set for each
  // MPI rank)
  ExecutionPolicy execPolicy =
      getExecutionPolicy(monteCarlo->processor_info->use_gpu);

  ParticleVaultContainer &my_particle_vault =
      *(monteCarlo->_particleVaultContainer);

  // Post Inital Receives for Particle Buffer
  monteCarlo->particle_buffer->Post_Receive_Particle_Buffer(
      my_particle_vault.getVaultSize());

  // Get Test For Done Method (Blocking or non-blocking
  MC_New_Test_Done_Method::Enum new_test_done_method =
      monteCarlo->particle_buffer->new_test_done_method;

  do 
  {
    int particle_count = 0; // Initialize count of num_particles processed

    while (!done) 
    {
      uint64_t fill_vault = 0;

      for (uint64_t processing_vault = 0;
           processing_vault < my_particle_vault.processingSize();
           processing_vault++) 
      {
        MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking_Kernel);
        uint64_t processed_vault =
            my_particle_vault.getFirstEmptyProcessedVault();

        ParticleVault *processingVault = my_particle_vault.getTaskProcessingVault(processing_vault);
        ParticleVault *processedVault = my_particle_vault.getTaskProcessedVault(processed_vault);

        int numParticles = processingVault->size();

        if (numParticles != 0) 
        {
          NVTX_Range trackingKernel(
              "cycleTracking_TrackingKernel"); // range ends at end of scope

          // The tracking kernel can run
          // * As a cuda kernel
          // * As an OpenMP 4.5 parallel loop on the GPU
          // * As an OpenMP 3.0 parallel loop on the CPU
          // * AS a single thread on the CPU.
          switch (execPolicy) {
          case gpuWithCUDA: {
#if defined(HAVE_CUDA)
            dim3 grid(1, 1, 1);
            dim3 block(1, 1, 1);
            int runKernel = ThreadBlockLayout(grid, block, numParticles);

            // Call Cycle Tracking Kernel
            if (runKernel)
              CycleTrackingKernel<<<grid, block>>>(
                  monteCarlo, numParticles, processingVault, processedVault);

            // Synchronize the stream so that memory is copied back before we
            // begin MPI section
            cudaPeekAtLastError();
            cudaDeviceSynchronize();
#endif
          } break;

          case gpuWithOpenMP: {
            int nthreads = 128;
            if (numParticles < 64 * 56)
              nthreads = 64;
            int nteams = (numParticles + nthreads - 1) / nthreads;
            nteams = nteams > 1 ? nteams : 1;
#ifdef HAVE_OPENMP_TARGET
#pragma omp target enter data map(to : monteCarlo [0:1])
#pragma omp target enter data map(to : processingVault [0:1])
#pragma omp target enter data map(to : processedVault [0:1])
#pragma omp target teams distribute parallel for num_teams(nteams)             \
    thread_limit(128)
#endif
            for (int particle_index = 0; particle_index < numParticles;
                 particle_index++) {
              CycleTrackingGuts(monteCarlo, particle_index, processingVault,
                                processedVault);
            }
#ifdef HAVE_OPENMP_TARGET
#pragma omp target exit data map(from : monteCarlo [0:1])
#pragma omp target exit data map(from : processingVault [0:1])
#pragma omp target exit data map(from : processedVault [0:1])
#endif
          } break;

          case cpu:
            #include "mc_omp_parallel_for_schedule_static.hh"
            for (int particle_index = 0; particle_index < numParticles;
                 particle_index++) 
            {
              // Tracking
              CycleTrackingGuts(monteCarlo, particle_index, processingVault,
                                processedVault);
            }
            break;
          default:
            qs_assert(false);
          } // end switch
        }

        particle_count += numParticles;

        MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking_Kernel);
        
        ARCANE_FATAL("youpi");

        MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking_MPI);

        // Next, communicate particles that have crossed onto
        // other MPI ranks.
        NVTX_Range cleanAndComm("cycleTracking_clean_and_comm");

        SendQueue &sendQueue = *(my_particle_vault.getSendQueue());
        monteCarlo->particle_buffer->Allocate_Send_Buffer(sendQueue);

        // Move particles from send queue to the send buffers
        for (int index = 0; index < sendQueue.size(); index++) 
        {
          sendQueueTuple &sendQueueT = sendQueue.getTuple(index);
          MC_Base_Particle mcb_particle;

          processingVault->getBaseParticleComm(mcb_particle,
                                               sendQueueT._particleIndex);

          int buffer = monteCarlo->particle_buffer->Choose_Buffer(sendQueueT._neighbor);
          monteCarlo->particle_buffer->Buffer_Particle(mcb_particle, buffer);
        }

        monteCarlo->particle_buffer->Send_Particle_Buffers(); // post MPI sends

        processingVault->clear(); // remove the invalid particles
        sendQueue.clear();

        // Move particles in "extra" vaults into the regular vaults.
        my_particle_vault.cleanExtraVaults();

        // receive any particles that have arrived from other ranks
        monteCarlo->particle_buffer->Receive_Particle_Buffers(fill_vault);

        MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking_MPI);

      } // for loop on vaults

      MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking_MPI);

      NVTX_Range collapseRange("cycleTracking_Collapse_ProcessingandProcessed");
      my_particle_vault.collapseProcessing();
      my_particle_vault.collapseProcessed();
      collapseRange.endRange();

      // Test for done - blocking on all MPI ranks
      NVTX_Range doneRange("cycleTracking_Test_Done_New");
      done = monteCarlo->particle_buffer->Test_Done_New(new_test_done_method);
      doneRange.endRange();

      MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking_MPI);

    } // while not done: Test_Done_New()

    // Everything should be done normally.
    done = monteCarlo->particle_buffer->Test_Done_New(
        MC_New_Test_Done_Method::Blocking);

  } while (!done);

  // Make sure to cancel all pending receive requests
  monteCarlo->particle_buffer->Cancel_Receive_Buffer_Requests();
  // Make sure Buffers Memory is Free
  monteCarlo->particle_buffer->Free_Buffers();
}

void QSModule::
trackingArc(MonteCarlo* monteCarlo)
{
  bool done = false;

  ParticleVaultContainer &my_particle_vault =
      *(monteCarlo->_particleVaultContainer);

  // Post Inital Receives for Particle Buffer
  monteCarlo->particle_buffer->Post_Receive_Particle_Buffer(
      my_particle_vault.getVaultSize());

  // Get Test For Done Method (Blocking or non-blocking
  MC_New_Test_Done_Method::Enum new_test_done_method =
      monteCarlo->particle_buffer->new_test_done_method;

  do 
  {
    int particle_count = 0; // Initialize count of num_particles processed

    while (!done) 
    {
      uint64_t fill_vault = 0;

      for (uint64_t processing_vault = 0;
           processing_vault < my_particle_vault.processingSize();
           processing_vault++) 
      {
        MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking_Kernel);
        uint64_t processed_vault =
            my_particle_vault.getFirstEmptyProcessedVault();

        ParticleVault *processingVault = my_particle_vault.getTaskProcessingVault(processing_vault);
        ParticleVault *processedVault = my_particle_vault.getTaskProcessedVault(processed_vault);

        int numParticles = processingVault->size();

        if (numParticles != 0) 
        {
          NVTX_Range trackingKernel(
              "cycleTracking_TrackingKernel"); // range ends at end of scope

          #include "mc_omp_parallel_for_schedule_static.hh"
          for (int particle_index = 0; particle_index < numParticles;
                particle_index++) 
          {
            //cout << particle_index;
            //flush(cout);
            // Tracking
            CycleTrackingGutsArc(monteCarlo, particle_index, processingVault,
                              processedVault);

          }
        }

        particle_count += numParticles;

        MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking_Kernel);

        ARCANE_FATAL("youpi");

        MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking_MPI);

        // Next, communicate particles that have crossed onto
        // other MPI ranks.
        NVTX_Range cleanAndComm("cycleTracking_clean_and_comm");

        SendQueue &sendQueue = *(my_particle_vault.getSendQueue());
        monteCarlo->particle_buffer->Allocate_Send_Buffer(sendQueue);

        // Move particles from send queue to the send buffers
        for (int index = 0; index < sendQueue.size(); index++) 
        {
          sendQueueTuple &sendQueueT = sendQueue.getTuple(index);
          MC_Base_Particle mcb_particle;

          processingVault->getBaseParticleComm(mcb_particle,
                                               sendQueueT._particleIndex);

          int buffer = monteCarlo->particle_buffer->Choose_Buffer(sendQueueT._neighbor);
          monteCarlo->particle_buffer->Buffer_Particle(mcb_particle, buffer);
        }

        monteCarlo->particle_buffer->Send_Particle_Buffers(); // post MPI sends

        processingVault->clear(); // remove the invalid particles
        sendQueue.clear();

        // Move particles in "extra" vaults into the regular vaults.
        my_particle_vault.cleanExtraVaults();

        // receive any particles that have arrived from other ranks
        monteCarlo->particle_buffer->Receive_Particle_Buffers(fill_vault);

        MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking_MPI);

      } // for loop on vaults

      MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking_MPI);

      NVTX_Range collapseRange("cycleTracking_Collapse_ProcessingandProcessed");
      my_particle_vault.collapseProcessing();
      my_particle_vault.collapseProcessed();
      collapseRange.endRange();

      // Test for done - blocking on all MPI ranks
      NVTX_Range doneRange("cycleTracking_Test_Done_New");
      done = monteCarlo->particle_buffer->Test_Done_New(new_test_done_method);
      doneRange.endRange();

      MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking_MPI);

    } // while not done: Test_Done_New()

    // Everything should be done normally.
    done = monteCarlo->particle_buffer->Test_Done_New(
        MC_New_Test_Done_Method::Blocking);

  } while (!done);

  // Make sure to cancel all pending receive requests
  monteCarlo->particle_buffer->Cancel_Receive_Buffer_Requests();
  // Make sure Buffers Memory is Free
  monteCarlo->particle_buffer->Free_Buffers();
}

void QSModule::
CycleTrackingGutsArc( MonteCarlo *monteCarlo, int particle_index, ParticleVault *processingVault, ParticleVault *processedVault )
{
    MC_Particle mc_particle;

    // Copy a single particle from the particle vault into mc_particle
    MC_Load_Particle(monteCarlo, mc_particle, processingVault, particle_index);

    // set the particle.task to the index of the processed vault the particle will census into.
    mc_particle.task = 0;//processed_vault;

    // loop over this particle until we cannot do anything more with it on this processor
    CycleTrackingFunctionArc( monteCarlo, mc_particle, particle_index, processingVault, processedVault );
    // if(particle_index < 11)
    // {
    // cout << "particle.identifier : " << mc_particle.identifier << endl;
    // cout << mc_particle.coordinate.x << " x " << mc_particle.coordinate.y << " x " << mc_particle.coordinate.z << endl;
    // }
    // else
    // {
    // exit(123);
    // }
    //Make sure this particle is marked as completed
    processingVault->invalidateParticle( particle_index );
}

void QSModule::
CycleTrackingFunctionArc( MonteCarlo *monteCarlo, MC_Particle &mc_particle, int particle_index, ParticleVault* processingVault, ParticleVault* processedVault)
{
    bool keepTrackingThisParticle = false;
    unsigned int tally_index = 0;
    unsigned int flux_tally_index = 0;
    unsigned int cell_tally_index = 0;
    //int DEBUG_compt = -1;
    do
    {
      //DEBUG_compt++;
        // Determine the outcome of a particle at the end of this segment such as:
        //
        //   (0) Undergo a collision within the current cell,
        //   (1) Cross a facet of the current cell,
        //   (2) Reach the end of the time step and enter census,
        //

        // Collision ou Census ou Facet crossing
        MC_Segment_Outcome_type::Enum segment_outcome = MC_Segment_OutcomeArc(monteCarlo, mc_particle, flux_tally_index);
        // if(particle_index == 11 && DEBUG_compt < 10)
        // {
        //   info() << "Passe ici !" << segment_outcome;
        // }

        m_numSegments = m_numSegments() + 1;

        mc_particle.num_segments += 1.;  /* Track the number of segments this particle has
                                            undergone this cycle on all processes. */
        switch (segment_outcome) {
        case MC_Segment_Outcome_type::Collision:
            {
            // The particle undergoes a collision event producing:
            //   (0) Other-than-one same-species secondary particle, or
            //   (1) Exactly one same-species secondary particle.
            if (CollisionEventArc(monteCarlo, mc_particle ) == MC_Collision_Event_Return::Continue_Tracking)
            {
                keepTrackingThisParticle = true;
            }
            else
            {
                keepTrackingThisParticle = false;
            }
            }
            break;
    
        case MC_Segment_Outcome_type::Facet_Crossing:
            {
                // The particle has reached a cell facet.
                MC_Tally_Event::Enum facet_crossing_type = MC_Facet_Crossing_EventArc(mc_particle, monteCarlo, particle_index, processingVault);
                // if(particle_index == 11 && DEBUG_compt < 10)
                // {
                //   info() << "   Passe ici !" << facet_crossing_type;
                // }
                if (facet_crossing_type == MC_Tally_Event::Facet_Crossing_Transit_Exit)
                {
                    keepTrackingThisParticle = true;  // Transit Event
                }
                else if (facet_crossing_type == MC_Tally_Event::Facet_Crossing_Escape)
                {
                    m_escape = m_escape() + 1;
                    mc_particle.last_event = MC_Tally_Event::Facet_Crossing_Escape;
                    mc_particle.species = -1;
                    keepTrackingThisParticle = false;
                }
                else if (facet_crossing_type == MC_Tally_Event::Facet_Crossing_Reflection)
                {
                    MCT_Reflect_ParticleArc(monteCarlo, mc_particle);
                    keepTrackingThisParticle = true;
                }
                else
                {
                    // Enters an adjacent cell in an off-processor domain.
                    //mc_particle.species = -1;
                    keepTrackingThisParticle = false;
                }
            }
            break;
    
        case MC_Segment_Outcome_type::Census:
            {
                // The particle has reached the end of the time step.
                processedVault->pushParticle(mc_particle);
                m_census = m_census() + 1;
                keepTrackingThisParticle = false;
                break;
            }
            
        default:
           qs_assert(false);
           break;  // should this be an error
        }
    
    } while ( keepTrackingThisParticle );
}

MC_Segment_Outcome_type::Enum QSModule::
MC_Segment_OutcomeArc(MonteCarlo* monteCarlo, MC_Particle &mc_particle, unsigned int &flux_tally_index)
{
    // initialize distances to large number
    int number_of_events = 3;
    double distance[3];
    distance[0] = distance[1] = distance[2] = 1e80; // +inf

    // Calculate the particle speed
    double particle_speed = mc_particle.Get_Velocity()->Length();

    // Force collision if a census event narrowly preempts a collision
    int force_collision = 0 ;
    if ( mc_particle.num_mean_free_paths < 0.0 )
    {
        force_collision = 1 ;

        if ( mc_particle.num_mean_free_paths > -900.0 )
        {
          printf(" MC_Segment_Outcome: mc_particle.num_mean_free_paths > -900.0 \n");
        }

        mc_particle.num_mean_free_paths = PhysicalConstants::_smallDouble;
    }

    // Randomly determine the distance to the next collision
    // based upon the composition of the current cell.
    double macroscopic_total_cross_section = weightedMacroscopicCrossSectionArc(monteCarlo, mc_particle.cellArc, mc_particle.energy_group);

    // Cache the cross section
    mc_particle.totalCrossSection = macroscopic_total_cross_section;
    if (macroscopic_total_cross_section == 0.0)
    {
        mc_particle.mean_free_path = PhysicalConstants::_hugeDouble;
    }
    else
    {
        mc_particle.mean_free_path = 1.0 / macroscopic_total_cross_section;
    }

    if ( mc_particle.num_mean_free_paths == 0.0)
    {
        // Sample the number of mean-free-paths remaining before
        // the next collision from an exponential distribution.
        double random_number = rngSample(&mc_particle.random_number_seed);

        mc_particle.num_mean_free_paths = -1.0*std::log(random_number);
    }

    // Calculate the distances to collision, nearest facet, and census.

    // Forced collisions do not need to move far.
    if (force_collision)
    {
        distance[MC_Segment_Outcome_type::Collision] = PhysicalConstants::_smallDouble;
    }
    else
    {
        distance[MC_Segment_Outcome_type::Collision] = mc_particle.num_mean_free_paths*mc_particle.mean_free_path;
    }

    // process census
    distance[MC_Segment_Outcome_type::Census] = particle_speed*mc_particle.time_to_census;


    //  DEBUG  Turn off threshold for now
    double distance_threshold = 10.0 * PhysicalConstants::_hugeDouble;
    // Get the current winning distance.
    double current_best_distance = PhysicalConstants::_hugeDouble;

    DirectionCosine *direction_cosine = mc_particle.Get_Direction_Cosine();

    bool new_segment =  (mc_particle.num_segments == 0 ||
                         mc_particle.last_event == MC_Tally_Event::Collision);

    // Calculate the minimum distance to each facet of the cell.
    MC_Nearest_Facet nearest_facet;
        nearest_facet = MCT_Nearest_FacetArc(&mc_particle, mc_particle.coordinate,
                                  direction_cosine, distance_threshold, current_best_distance, new_segment, monteCarlo);

    Face faceOfFacet = mc_particle.cellArc.face(nearest_facet.facet / 4);

    mc_particle.normal_dot = nearest_facet.dot_product;

    distance[MC_Segment_Outcome_type::Facet_Crossing] = nearest_facet.distance_to_facet;
    //info() << "Distance : " << nearest_facet.distance_to_facet << " " << nearest_facet.facet;
    //ARCANE_FATAL("aaa");

    // Get out of here if the tracker failed to bound this particle's volume.
    if (mc_particle.last_event == MC_Tally_Event::Facet_Crossing_Tracking_Error)
    {
        return MC_Segment_Outcome_type::Facet_Crossing;
    }

    // Calculate the minimum distance to the selected events.

    // Force a collision (if required).
    if ( force_collision == 1 )
    {
        distance[MC_Segment_Outcome_type::Facet_Crossing] = PhysicalConstants::_hugeDouble;
        distance[MC_Segment_Outcome_type::Census]         = PhysicalConstants::_hugeDouble;
        distance[MC_Segment_Outcome_type::Collision]      = PhysicalConstants::_tinyDouble ;
    }

    // we choose our segment outcome here
    MC_Segment_Outcome_type::Enum segment_outcome =
        (MC_Segment_Outcome_type::Enum) MC_Find_Min(distance, number_of_events);

    if (distance[segment_outcome] < 0)
    {
        MC_Fatal_Jump( "Negative distances to events are NOT permitted!\n"
                       "identifier              = %" PRIu64 "\n"
                       "(Collision              = %g,\n"
                       " Facet Crossing         = %g,\n"
                       " Census                 = %g,\n",
                       mc_particle.identifier,
                       distance[MC_Segment_Outcome_type::Collision],
                       distance[MC_Segment_Outcome_type::Facet_Crossing],
                       distance[MC_Segment_Outcome_type::Census]);
    }
    mc_particle.segment_path_length = distance[segment_outcome];

    mc_particle.num_mean_free_paths -= mc_particle.segment_path_length / mc_particle.mean_free_path;

    // Before using segment_outcome as an index, verify it is valid
    if (segment_outcome < 0 || segment_outcome >= MC_Segment_Outcome_type::Max_Number)
    {
        MC_Fatal_Jump( "segment_outcome '%d' is invalid\n", (int)segment_outcome );
    }

    MC_Tally_Event::Enum SegmentOutcome_to_LastEvent[MC_Segment_Outcome_type::Max_Number] =
    {
        MC_Tally_Event::Collision,
        MC_Tally_Event::Facet_Crossing_Transit_Exit,
        MC_Tally_Event::Census,
    };

    mc_particle.last_event = SegmentOutcome_to_LastEvent[segment_outcome];

    // Set the segment path length to be the minimum of
    //   (i)   the distance to collision in the cell, or
    //   (ii)  the minimum distance to a facet of the cell, or
    //   (iii) the distance to census at the end of the time step
    if (segment_outcome == MC_Segment_Outcome_type::Collision)
    {
        mc_particle.num_mean_free_paths = 0.0;
    }
    else if (segment_outcome == MC_Segment_Outcome_type::Facet_Crossing)
    {
        mc_particle.facet = nearest_facet.facet;
        mc_particle.faceArc = faceOfFacet;
    }
    else if (segment_outcome == MC_Segment_Outcome_type::Census)
    {
        mc_particle.time_to_census = MC_MIN(mc_particle.time_to_census, 0.0);
    }

    // If collision was forced, set mc_particle.num_mean_free_paths = 0
    // so that a new value is randomly selected on next pass.
    if (force_collision == 1) { mc_particle.num_mean_free_paths = 0.0; }

    // Do not perform any tallies if the segment path length is zero.
    //   This only introduces roundoff errors.
    if (mc_particle.segment_path_length == 0.0)
    {
        return segment_outcome;
    }

    // Move particle to end of segment, accounting for some physics processes along the segment.

    // Project the particle trajectory along the segment path length.
    mc_particle.Move_Particle(mc_particle.direction_cosine, mc_particle.segment_path_length);

    double segment_path_time = (mc_particle.segment_path_length/particle_speed);

    // Decrement the time to census and increment age.
    mc_particle.time_to_census -= segment_path_time;
    mc_particle.age += segment_path_time;

      //info() << "Drap" << mc_particle.age;
    // Ensure mc_particle.time_to_census is non-negative.
    if (mc_particle.time_to_census < 0.0)
    {
        mc_particle.time_to_census = 0.0;
    }

    // Accumulate the particle's contribution to the scalar flux.
    m_scalarFluxTally[mc_particle.cellArc][mc_particle.energy_group] += mc_particle.segment_path_length * mc_particle.weight;

    return segment_outcome;
}

double QSModule::
weightedMacroscopicCrossSectionArc(MonteCarlo* monteCarlo, Cell &cell, int energyGroup)
{
   double precomputedCrossSection = m_total[cell][energyGroup];

   if (precomputedCrossSection > 0.0)
      return precomputedCrossSection;
   
   int globalMatIndex = m_material[cell];
   int nIsotopes = (int)monteCarlo->_materialDatabase->_mat[globalMatIndex]._iso.size();
   double sum = 0.0;
   for (int isoIndex = 0; isoIndex < nIsotopes; isoIndex++)
   {
      sum += macroscopicCrossSectionArc(monteCarlo, -1, cell, isoIndex, energyGroup);
   }

   ATOMIC_WRITE( m_total[cell][energyGroup], sum );

   return sum;
}

//----------------------------------------------------------------------------------------------------------------------
//  Routine MacroscopicCrossSection calculates the number-density-weighted macroscopic cross
//  section of a cell.
//
//  A reactionIndex of -1 means total cross section.
//----------------------------------------------------------------------------------------------------------------------

double QSModule::
macroscopicCrossSectionArc(MonteCarlo* monteCarlo, int reactionIndex, Cell& cell, int isoIndex, int energyGroup)
{
   // Initialize various data items.
   int globalMatIndex = m_material[cell];

   double atomFraction = monteCarlo->_materialDatabase->_mat[globalMatIndex]._iso[isoIndex]._atomFraction;

   double microscopicCrossSection = 0.0;
   // The cell number density is the fraction of the atoms in cell
   // volume of this isotope.  We set this (elsewhere) to 1/nIsotopes.
   // This is a statement that we treat materials as if all of their
   // isotopes are present in equal amounts
   double cellNumberDensity = m_cellNumberDensity[cell];

   int isotopeGid = monteCarlo->_materialDatabase->_mat[globalMatIndex]._iso[isoIndex]._gid;
   if ( atomFraction == 0.0 || cellNumberDensity == 0.0) { return 1e-20; }

   if (reactionIndex < 0)
   {
      // Return total cross section
      microscopicCrossSection = monteCarlo->_nuclearData->getTotalCrossSection(isotopeGid, energyGroup);
   }
   else
   {
      // Return the reaction cross section
      microscopicCrossSection = monteCarlo->_nuclearData->getReactionCrossSection((unsigned int)reactionIndex,
                isotopeGid, energyGroup);
   }

   return atomFraction * cellNumberDensity * microscopicCrossSection;

}

MC_Nearest_Facet QSModule::
MCT_Nearest_FacetArc( MC_Particle *mc_particle,
                      MC_Vector &coordinate,
                      const DirectionCosine *direction_cosine,
                      double distance_threshold,
                      double current_best_distance,
                      bool new_segment,
                      MonteCarlo* monteCarlo )
{

    MC_Nearest_Facet nearest_facet =
       MCT_Nearest_Facet_3D_GArc(mc_particle, coordinate, direction_cosine);

    if (nearest_facet.distance_to_facet < 0) {
      nearest_facet.distance_to_facet = 0; 
    }

    if (nearest_facet.distance_to_facet >= PhysicalConstants::_hugeDouble)
    {
       qs_assert(false);
    }

    return nearest_facet;
}  // End MCT_Nearest_Facet

MC_Nearest_Facet QSModule::
MCT_Nearest_Facet_3D_GArc(MC_Particle *mc_particle,
                          MC_Vector &coordinate,
                          const DirectionCosine *direction_cosine)
{
  Cell cell = mc_particle->cellArc;
  MC_Vector *facet_coords[3];
  int iteration = 0;
  double move_factor = 0.5 * PhysicalConstants::_smallDouble;

  // Initialize some data for the unstructured, hexahedral mesh.
  int num_facets_per_cell = 24;

  while (true) // will break out when distance is found
  {
      // Determine the distance to each facet of the cell.
      // (1e-8 * Radius)^2
      double plane_tolerance = 1e-16*(coordinate.x*coordinate.x +
                                      coordinate.y*coordinate.y +
                                      coordinate.z*coordinate.z);

      MC_Distance_To_Facet distance_to_facet[24];

      int facet_index = -1;
      ENUMERATE_FACE(iface, cell.faces())
      {
        Face face = *iface;
        

        for (int i = 0; i < 4; i++)
        {
          facet_index++;
          if((int)facet_index/4 != iface.index())
          {
            ARCANE_FATAL("Erreur facet index");
          }
          int first_pos_node = (ordre_qs[iface.index()] ? ((i == 3) ? 0 : i+1) : i);
          int second_pos_node = (ordre_qs[iface.index()] ? i : ((i == 3) ? 0 : i+1));

          Node first_node = face.node(first_pos_node);
          Node second_node = face.node(second_pos_node);

          MC_Vector point0(m_coordCm[first_node][0], m_coordCm[first_node][1], m_coordCm[first_node][2]);
          MC_Vector point1(m_coordCm[second_node][0], m_coordCm[second_node][1], m_coordCm[second_node][2]);
          MC_Vector point2(m_coordMidCm[iface][0], m_coordMidCm[iface][1], m_coordMidCm[iface][2]);

          distance_to_facet[facet_index].distance = PhysicalConstants::_hugeDouble;

          MC_General_Plane plane(point0, point1, point2);

          double facet_normal_dot_direction_cosine =
              (plane.A * direction_cosine->alpha +
              plane.B * direction_cosine->beta +
              plane.C * direction_cosine->gamma);

          // info() << " Facet : " << point0.x << "x" << point0.y << "x" << point0.z;
          // info()                << point1.x << "x" << point1.y << "x" << point1.z;
          // info()                << point2.x << "x" << point2.y << "x" << point2.z;
          // info() << " Coord : " << m_coord[first_node][0] << "x" << m_coord[first_node][1] << "x" << m_coord[first_node][2];
          // info()                << m_coord[second_node][0] << "x" << m_coord[second_node][1] << "x" << m_coord[second_node][2];
          // info()                << m_coordMid[iface][0] << "x" << m_coordMid[iface][1] << "x" << m_coordMid[iface][2];
          //info() << " CooCm : " << m_coordCm[first_node][0] << "x" << m_coordCm[first_node][1] << "x" << m_coordCm[first_node][2];
          //info()                << m_coordCm[second_node][0] << "x" << m_coordCm[second_node][1] << "x" << m_coordCm[second_node][2];
          //info()                << m_coordMidCm[iface][0] << "x" << m_coordMidCm[iface][1] << "x" << m_coordMidCm[iface][2];
          //info() << " Center : " << m_coordCenter[cell][0] << "x" << m_coordCenter[cell][1] << "x" << m_coordCenter[cell][2];
          //info() << " Plane : " << plane.A << "x" << plane.B << "x" << plane.C << " " << plane.D;
          //info() << "facet_normal_dot_direction_cosine : " << facet_normal_dot_direction_cosine;
          // info() << " Facet : " << m_coordCm[first_node][0] << "x" << m_coordCm[first_node][1] << "x" << m_coordCm[first_node][2];
          // info()                << m_coordCm[second_node][0] << "x" << m_coordCm[second_node][1] << "x" << m_coordCm[second_node][2];
          // info()                << m_coordMidCm[iface][0] << "x" << m_coordMidCm[iface][1] << "x" << m_coordMidCm[iface][2];
          // info();
          // info();
          //info() << "Coord : " << coordinate.x << "x" << coordinate.y << "x" << coordinate.z;
          

          // Consider only those facets whose outer normals have
          // a positive dot product with the direction cosine.
          // I.e. the particle is LEAVING the cell.
          if (facet_normal_dot_direction_cosine <= 0.0) { continue; }

          double t = MCT_Nearest_Facet_3D_G_Distance_To_Segment(
              plane_tolerance,
              facet_normal_dot_direction_cosine, plane.A, plane.B, plane.C, plane.D,
              point0, point1, point2,
              coordinate, direction_cosine, false);

          distance_to_facet[facet_index].distance = t;
          //info() << "distance_to_facet : facet_index :" << facet_index << " distance : " << t;
          //info() << "mcp :" << mc_particle->identifier;
          //info() << "direction_cosine : " << direction_cosine->alpha << " x "<< direction_cosine->beta << " x "<< direction_cosine->gamma;
          //info() << "\n";

        }
      }
      //ARCANE_FATAL("aaa");

      int retry = 0;

      MC_Nearest_Facet nearest_facet = MCT_Nearest_Facet_Find_NearestArc(
        mc_particle, coordinate,
        iteration, move_factor, num_facets_per_cell,
//to-do       monteCarlo->distance_to_facet->task[my_task_num].facet,
        distance_to_facet,
        retry);


      if (! retry) return nearest_facet;
  } // while (true)
}  // End MCT_Nearest_Facet_3D_G

double QSModule::
MCT_Nearest_Facet_3D_G_Distance_To_Segment(double plane_tolerance,
                                                     double facet_normal_dot_direction_cosine,//=
                                                     double A, double B, double C, double D,//=
                                                     const MC_Vector &facet_coords0,//=
                                                     const MC_Vector &facet_coords1,//=
                                                     const MC_Vector &facet_coords2,//=
                                                     const MC_Vector &coordinate,
                                                     const DirectionCosine *direction_cosine,//=
                                                     bool allow_enter) //=
   {
    double boundingBox_tolerance = 1e-9;
    double numerator = -1.0*(A * coordinate.x +
                             B * coordinate.y +
                             C * coordinate.z +
                             D);

    /* Plane equation: numerator = -P(x,y,z) = -(Ax + By + Cz + D)
       if: numerator < -1e-8*length(x,y,z)   too negative!
       if: numerator < 0 && numerator^2 > ( 1e-8*length(x,y,z) )^2   too negative!
       reverse inequality since squaring function is decreasing for negative inputs.
       If numerator is just SLIGHTLY negative, then the particle is just outside of the face */

    // Filter out too negative distances
    if (!allow_enter && numerator < 0.0 && numerator * numerator > plane_tolerance) {
        return PhysicalConstants::_hugeDouble; }

    // we have to restrict the solution to within the triangular face
    double distance = numerator / facet_normal_dot_direction_cosine;

    // see if the intersection point of the ray and the plane is within the triangular facet
    MC_Vector intersection_pt;
    intersection_pt.x = coordinate.x + distance * direction_cosine->alpha;
    intersection_pt.y = coordinate.y + distance * direction_cosine->beta;
    intersection_pt.z = coordinate.z + distance * direction_cosine->gamma;

    // if the point is completely below the triangle, it is not in the triangle
#define IF_POINT_BELOW_CONTINUE(axis)                                    \
        if ( facet_coords0.axis > intersection_pt.axis + boundingBox_tolerance&& \
             facet_coords1.axis > intersection_pt.axis + boundingBox_tolerance && \
             facet_coords2.axis > intersection_pt.axis + boundingBox_tolerance ) { return PhysicalConstants::_hugeDouble; }

#define IF_POINT_ABOVE_CONTINUE(axis)                                    \
        if ( facet_coords0.axis < intersection_pt.axis - boundingBox_tolerance && \
             facet_coords1.axis < intersection_pt.axis - boundingBox_tolerance && \
             facet_coords2.axis < intersection_pt.axis - boundingBox_tolerance ) { return PhysicalConstants::_hugeDouble; }

    // Is the intersection point inside the triangular facet?  Project to 2D and see.

    // A^2 + B^2 + C^2 = 1, so max(|A|,|B|,|C|) >= 1/sqrt(3) = 0.577
    // (all coefficients can't be small)
    double cross0 = 0, cross1 = 0, cross2 = 0;
    if ( C < -0.5 || C > 0.5 )
    {
       IF_POINT_BELOW_CONTINUE(x);
       IF_POINT_ABOVE_CONTINUE(x);
       IF_POINT_BELOW_CONTINUE(y);
       IF_POINT_ABOVE_CONTINUE(y);

#define AB_CROSS_AC(ax,ay,bx,by,cx,cy) ( (bx-ax)*(cy-ay) - (by-ay)*(cx-ax) )

       cross1 = AB_CROSS_AC(facet_coords0.x, facet_coords0.y,
                            facet_coords1.x, facet_coords1.y,
                            intersection_pt.x,  intersection_pt.y);
       cross2 = AB_CROSS_AC(facet_coords1.x, facet_coords1.y,
                            facet_coords2.x, facet_coords2.y,
                            intersection_pt.x,  intersection_pt.y);
       cross0 = AB_CROSS_AC(facet_coords2.x, facet_coords2.y,
                            facet_coords0.x, facet_coords0.y,
                            intersection_pt.x,  intersection_pt.y);

    }
    else if ( B < -0.5 || B > 0.5 )
    {
       IF_POINT_BELOW_CONTINUE(x);
       IF_POINT_ABOVE_CONTINUE(x);
       IF_POINT_BELOW_CONTINUE(z);
       IF_POINT_ABOVE_CONTINUE(z);

       cross1 = AB_CROSS_AC(facet_coords0.z, facet_coords0.x,
                            facet_coords1.z, facet_coords1.x,
                            intersection_pt.z,  intersection_pt.x);
       cross2 = AB_CROSS_AC(facet_coords1.z, facet_coords1.x,
                            facet_coords2.z, facet_coords2.x,
                            intersection_pt.z,  intersection_pt.x);
       cross0 = AB_CROSS_AC(facet_coords2.z, facet_coords2.x,
                            facet_coords0.z, facet_coords0.x,
                            intersection_pt.z,  intersection_pt.x);

    }
    else if ( A < -0.5 || A > 0.5 )
    {
       IF_POINT_BELOW_CONTINUE(z);
       IF_POINT_ABOVE_CONTINUE(z);
       IF_POINT_BELOW_CONTINUE(y);
       IF_POINT_ABOVE_CONTINUE(y);

       cross1 = AB_CROSS_AC(facet_coords0.y, facet_coords0.z,
                            facet_coords1.y, facet_coords1.z,
                            intersection_pt.y,  intersection_pt.z);
       cross2 = AB_CROSS_AC(facet_coords1.y, facet_coords1.z,
                            facet_coords2.y, facet_coords2.z,
                            intersection_pt.y,  intersection_pt.z);
       cross0 = AB_CROSS_AC(facet_coords2.y, facet_coords2.z,
                            facet_coords0.y, facet_coords0.z,
                            intersection_pt.y,  intersection_pt.z);
    }

    double cross_tol = 1e-9 * MC_FABS(cross0 + cross1 + cross2);  // cross product tolerance

    if ( (cross0 > -cross_tol && cross1 > -cross_tol && cross2 > -cross_tol) ||
         (cross0 <  cross_tol && cross1 <  cross_tol && cross2 <  cross_tol) )
    {
        return distance;
    }
    return PhysicalConstants::_hugeDouble;
   }

MC_Nearest_Facet QSModule::
MCT_Nearest_Facet_Find_NearestArc(MC_Particle *mc_particle,
                                  MC_Vector &coordinate,
                                  int &iteration, // input/output
                                  double &move_factor, // input/output
                                  int num_facets_per_cell,
                                  MC_Distance_To_Facet *distance_to_facet,
                                  int &retry /* output */ )
{
  MC_Nearest_Facet nearest_facet = MCT_Nearest_Facet_Find_Nearest(num_facets_per_cell, distance_to_facet);

  const int max_allowed_segments = 10000000;

  retry = 0;

  if ( mc_particle )
  {
      if ( (nearest_facet.distance_to_facet == PhysicalConstants::_hugeDouble && move_factor > 0) ||
          ( mc_particle->num_segments > max_allowed_segments && nearest_facet.distance_to_facet <= 0.0 ) )
      {
        // Could not find a solution, so move the particle towards the center of the cell
        // and try again.
        MCT_Nearest_Facet_3D_G_Move_ParticleArc(mc_particle->cellArc, coordinate, move_factor);

        iteration++;
        move_factor *= 2.0;

        if ( move_factor > 1.0e-2 )
            move_factor = 1.0e-2;

        int max_iterations = 10000;

        if ( iteration == max_iterations )
        {
          //info() << (nearest_facet.distance_to_facet == PhysicalConstants::_hugeDouble) << (move_factor > 0) << 
          //(mc_particle->num_segments > max_allowed_segments) << (nearest_facet.distance_to_facet <= 0.0);

            qs_assert(false); // If we start hitting this assertion we can
            // come up with a better mitigation plan. - dfr
            retry = 0;

        }
        else
            retry = 1;

        // Allow the distance to the current facet
        //location->facet = -1;

      }
  }
  return nearest_facet;
}

MC_Nearest_Facet QSModule::
MCT_Nearest_Facet_Find_Nearest( int num_facets_per_cell,
                                MC_Distance_To_Facet *distance_to_facet)
   {
      MC_Nearest_Facet nearest_facet;

      // largest negative distance (smallest magnitude, but negative)
      MC_Nearest_Facet nearest_negative_facet;
      nearest_negative_facet.distance_to_facet = -PhysicalConstants::_hugeDouble;

      // Determine the facet that is closest to the specified coordinates.
      for (int facet_index = 0; facet_index < num_facets_per_cell; facet_index++)
      {
         if ( distance_to_facet[facet_index].distance > 0.0 )
         {
            if ( distance_to_facet[facet_index].distance <= nearest_facet.distance_to_facet )
            {
               nearest_facet.distance_to_facet = distance_to_facet[facet_index].distance;
               nearest_facet.facet             = facet_index;
            }
         }
         else // zero or negative distance
         {
            if ( distance_to_facet[facet_index].distance > nearest_negative_facet.distance_to_facet )
            {
               // smallest in magnitude, but negative
               nearest_negative_facet.distance_to_facet = distance_to_facet[facet_index].distance;
               nearest_negative_facet.facet             = facet_index;
            }
         }
      }


      if ( nearest_facet.distance_to_facet == PhysicalConstants::_hugeDouble )
      {
         if ( nearest_negative_facet.distance_to_facet != -PhysicalConstants::_hugeDouble )
         {
            // no positive solution, so allow a negative solution, that had really small magnitude.
            nearest_facet.distance_to_facet = nearest_negative_facet.distance_to_facet;
            nearest_facet.facet             = nearest_negative_facet.facet;
         }
      }

      return nearest_facet;
   }

void QSModule::
MCT_Nearest_Facet_3D_G_Move_ParticleArc( Cell& cell,
                                      MC_Vector &coordinate, // input/output: move this coordinate
                                      double move_factor)      // input: multiplication factor for move
   {
      MC_Vector move_to(m_coordCenter[cell][0], m_coordCenter[cell][1], m_coordCenter[cell][2]);

      coordinate.x += move_factor * ( move_to.x - coordinate.x );
      coordinate.y += move_factor * ( move_to.y - coordinate.y );
      coordinate.z += move_factor * ( move_to.z - coordinate.z );
   }

unsigned int QSModule::
MC_Find_Min(const double *array, int num_elements)
{
    double min = array[0];
    int    min_index = 0;

    for (int element_index = 1; element_index < num_elements; ++element_index)
    {
        if ( array[element_index] < min )
        {
            min = array[element_index];
            min_index = element_index;
        }
    }

    return min_index;
}

bool QSModule::
CollisionEventArc(MonteCarlo* monteCarlo, MC_Particle &mc_particle)
{
  Cell cell = mc_particle.cellArc;
   int globalMatIndex = m_material[cell];

   //------------------------------------------------------------------------------------------------------------------
   //    Pick the isotope and reaction.
   //------------------------------------------------------------------------------------------------------------------
   double randomNumber = rngSample(&mc_particle.random_number_seed);
   double totalCrossSection = mc_particle.totalCrossSection;
   double currentCrossSection = totalCrossSection * randomNumber;
   int selectedIso = -1;
   int selectedUniqueNumber = -1;
   int selectedReact = -1;
   int numIsos = (int)monteCarlo->_materialDatabase->_mat[globalMatIndex]._iso.size();
   
   for (int isoIndex = 0; isoIndex < numIsos && currentCrossSection >= 0; isoIndex++)
   {
      int uniqueNumber = monteCarlo->_materialDatabase->_mat[globalMatIndex]._iso[isoIndex]._gid;
      int numReacts = monteCarlo->_nuclearData->getNumberReactions(uniqueNumber);
      for (int reactIndex = 0; reactIndex < numReacts; reactIndex++)
      {
         currentCrossSection -= macroscopicCrossSectionArc(monteCarlo, reactIndex, mc_particle.cellArc,
                   isoIndex, mc_particle.energy_group);
         if (currentCrossSection < 0)
         {
            selectedIso = isoIndex;
            selectedUniqueNumber = uniqueNumber;
            selectedReact = reactIndex;
            break;
         }
      }
   }
   qs_assert(selectedIso != -1);

   //------------------------------------------------------------------------------------------------------------------
   //    Do the collision.
   //------------------------------------------------------------------------------------------------------------------
   double energyOut[MAX_PRODUCTION_SIZE];
   double angleOut[MAX_PRODUCTION_SIZE];
   int nOut = 0;
   double mat_mass = monteCarlo->_materialDatabase->_mat[globalMatIndex]._mass;

   monteCarlo->_nuclearData->_isotopes[selectedUniqueNumber]._species[0]._reactions[selectedReact].sampleCollision(
      mc_particle.kinetic_energy, mat_mass, &energyOut[0], &angleOut[0], nOut, &(mc_particle.random_number_seed), MAX_PRODUCTION_SIZE );

   //--------------------------------------------------------------------------------------------------------------
   //  Post-Collision Phase 1:
   //    Tally the collision
   //--------------------------------------------------------------------------------------------------------------

   // Set the reaction for this particle.
   m_collision = m_collision() + 1;
   NuclearDataReaction::Enum reactionType = monteCarlo->_nuclearData->_isotopes[selectedUniqueNumber]._species[0].\
           _reactions[selectedReact]._reactionType;

   switch (reactionType)
   {
      case NuclearDataReaction::Scatter:
         m_scatter = m_scatter() + 1;
         break;
      case NuclearDataReaction::Absorption:
         m_absorb = m_absorb() + 1;
         break;
      case NuclearDataReaction::Fission:
         m_fission = m_fission() + 1;
         m_produce = m_produce() + nOut;
         break;
      case NuclearDataReaction::Undefined:
         printf("reactionType invalid\n");
         qs_assert(false);
   }

   if( nOut == 0 ) return false;

   for (int secondaryIndex = 1; secondaryIndex < nOut; secondaryIndex++)
   {
        // Newly created particles start as copies of their parent
        MC_Particle secondaryParticle = mc_particle;
        secondaryParticle.random_number_seed = rngSpawn_Random_Number_Seed(&mc_particle.random_number_seed);
        secondaryParticle.identifier = secondaryParticle.random_number_seed;
        updateTrajectory( energyOut[secondaryIndex], angleOut[secondaryIndex], secondaryParticle );
        monteCarlo->_particleVaultContainer->addExtraParticle(secondaryParticle);
   }

   updateTrajectory( energyOut[0], angleOut[0], mc_particle);

   // If a fission reaction produces secondary particles we also add the original
   // particle to the "extras" that we will handle later.  This avoids the 
   // possibility of a particle doing multiple fission reactions in a single
   // kernel invocation and overflowing the extra storage with secondary particles.
   if ( nOut > 1 ) 
       monteCarlo->_particleVaultContainer->addExtraParticle(mc_particle);

   //If we are still tracking this particle the update its energy group
   mc_particle.energy_group = monteCarlo->_nuclearData->getEnergyGroup(mc_particle.kinetic_energy);

   return nOut == 1;
}

void QSModule::
updateTrajectory( double energy, double angle, MC_Particle& particle )
{
    particle.kinetic_energy = energy;
    double cosTheta = angle;
    double randomNumber = rngSample(&particle.random_number_seed);
    double phi = 2 * 3.14159265 * randomNumber;
    double sinPhi = sin(phi);
    double cosPhi = cos(phi);
    double sinTheta = sqrt((1.0 - (cosTheta*cosTheta)));
    particle.direction_cosine.Rotate3DVector(sinTheta, cosTheta, sinPhi, cosPhi);
    double speed = (PhysicalConstants::_speedOfLight *
            sqrt((1.0 - ((PhysicalConstants::_neutronRestMassEnergy *
            PhysicalConstants::_neutronRestMassEnergy) /
            ((energy + PhysicalConstants::_neutronRestMassEnergy) *
            (energy + PhysicalConstants::_neutronRestMassEnergy))))));
    particle.velocity.x = speed * particle.direction_cosine.alpha;
    particle.velocity.y = speed * particle.direction_cosine.beta;
    particle.velocity.z = speed * particle.direction_cosine.gamma;
    randomNumber = rngSample(&particle.random_number_seed);
    particle.num_mean_free_paths = -1.0*std::log(randomNumber);
}

MC_Tally_Event::Enum QSModule::
MC_Facet_Crossing_EventArc(MC_Particle &mc_particle, MonteCarlo* monteCarlo, int particle_index, ParticleVault* processingVault)
{
     //Subfacet_Adjacency &facet_adjacency = MCT_Adjacent_Facet(location, mc_particle, monteCarlo);
    Face face = mc_particle.faceArc;

    if ( m_boundaryCond[face] == MC_Subfacet_Adjacency_Event::Transit_On_Processor )
    {
        // The particle will enter into an adjacent cell.
        //mc_particle.domain     = facet_adjacency.adjacent.domain;
        Cell cell = face.frontCell();

        if(cell == mc_particle.cellArc)
        {
          cell = face.backCell();
        }
        if(cell == mc_particle.cellArc)
        {
          if(face.frontCell() == face.backCell())
          {
            ARCANE_FATAL("hhhhhhhhhhhaaaaaaaaaaaaaa");
          }
          ARCANE_FATAL("Erreur deplace");
        }
        mc_particle.cellArc = cell;
        mc_particle.cell       = cell.localId();
        mc_particle.facet      = (mc_particle.facet < 12 ? mc_particle.facet + 12 : mc_particle.facet - 12);
        mc_particle.last_event = MC_Tally_Event::Facet_Crossing_Transit_Exit;
    }
    else if ( m_boundaryCond[face] == MC_Subfacet_Adjacency_Event::Boundary_Escape )
    {
        // The particle will escape across the system boundary.
        mc_particle.last_event = MC_Tally_Event::Facet_Crossing_Escape;
    }
    else if ( m_boundaryCond[face] == MC_Subfacet_Adjacency_Event::Boundary_Reflection )
    {
        // The particle will reflect off of the system boundary.
        mc_particle.last_event = MC_Tally_Event::Facet_Crossing_Reflection;
    }
    else if ( m_boundaryCond[face] == MC_Subfacet_Adjacency_Event::Transit_Off_Processor )
    {
        // The particle will enter into an adjacent cell on a spatial neighbor.
        // The neighboring domain is on another processor. Set domain local domain on neighbor proc
        
        //mc_particle.domain     = facet_adjacency.adjacent.domain;
        ARCANE_FATAL("TODO non supporté");
        int id_face = mc_particle.facet / 4;
        Cell cell;
        if(id_face < 3)
        {
          cell = face.backCell();
        }
        else
        {
          cell = face.frontCell();
        }
        mc_particle.cellArc = cell;
        mc_particle.facet      = (mc_particle.facet < 12 ? mc_particle.facet + 12 : mc_particle.facet - 12);
        mc_particle.last_event = MC_Tally_Event::Facet_Crossing_Communication;

        // Select particle buffer
        // TODO Gérer transfert particles.
        //int neighbor_rank = monteCarlo->domain[facet_adjacency.current.domain].mesh._nbrRank[facet_adjacency.neighbor_index];

        processingVault->putParticle( mc_particle, particle_index );

        //Push neighbor rank and mc_particle onto the send queue
        //monteCarlo->_particleVaultContainer->getSendQueue()->push( neighbor_rank, particle_index );

    }

    return mc_particle.last_event;
}

void QSModule::
MCT_Reflect_ParticleArc(MonteCarlo *monteCarlo, MC_Particle &particle)
{
    DirectionCosine *direction_cosine = particle.Get_Direction_Cosine();

    int truc = particle.facet % 4;

    Cell cell = particle.cellArc;
    Face face = particle.faceArc;

    int first_pos_node = (ordre_qs[m_indexArc[face]] ? ((truc == 3) ? 0 : truc+1) : truc);
    int second_pos_node = (ordre_qs[m_indexArc[face]] ? truc : ((truc == 3) ? 0 : truc+1));

    Node first_node = face.node(first_pos_node);
    Node second_node = face.node(second_pos_node);


    MC_Vector point0(m_coordCm[first_node][0], m_coordCm[first_node][1], m_coordCm[first_node][2]);
    MC_Vector point1(m_coordCm[second_node][0], m_coordCm[second_node][1], m_coordCm[second_node][2]);
    MC_Vector point2(m_coordMidCm[face][0], m_coordMidCm[face][1], m_coordMidCm[face][2]);


    MC_General_Plane plane(point0, point1, point2);

    MC_Vector facet_normal(plane.A, plane.B, plane.C);


    double dot = 2.0*( direction_cosine->alpha * facet_normal.x +
                       direction_cosine->beta  * facet_normal.y +
                       direction_cosine->gamma * facet_normal.z );

    if ( dot > 0 ) // do not reflect a particle that is ALREADY pointing inward
    {
        // reflect the particle
        direction_cosine->alpha -= dot * facet_normal.x;
        direction_cosine->beta  -= dot * facet_normal.y;
        direction_cosine->gamma -= dot * facet_normal.z;
    }

    // Calculate the reflected, velocity components.
    double particle_speed = particle.velocity.Length();
    particle.velocity.x = particle_speed * particle.direction_cosine.alpha;
    particle.velocity.y = particle_speed * particle.direction_cosine.beta;
    particle.velocity.z = particle_speed * particle.direction_cosine.gamma;
}
