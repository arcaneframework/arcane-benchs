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

#include "arcane/IParallelMng.h"
#include "arcane/IMesh.h"

using namespace Arcane;
using namespace std;

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

  monteCarlo->_tallies->CycleInitialize(monteCarlo);

  monteCarlo->_particleVaultContainer->swapProcessingProcessedVaults();
  monteCarlo->_particleVaultContainer->collapseProcessed();
  monteCarlo->_particleVaultContainer->collapseProcessing();

  // Nombre de particles dans Processing.
  monteCarlo->_tallies->_balanceTask[0]._start =
      monteCarlo->_particleVaultContainer->sizeProcessing();

  monteCarlo->particle_buffer->Initialize();

  // Création des particules (pas compris la répartition du travail entre proc).
  MC_SourceNow(monteCarlo);

  // Réduction ou augmentation du nombre de particules.
  PopulationControl(monteCarlo, loadBalance); // controls particle population

  // Roulette sur les particules avec faible poids.
  RouletteLowWeightParticles(monteCarlo); // Delete particles with low statistical weight

  MC_FASTTIMER_STOP(MC_Fast_Timer::cycleInit);
}

void QSModule::
cycleTracking()
{
  info() << "Module Quicksilver cycleTracking";

  MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking);

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
    m_coord.resize(4);
  #else
    m_coord.resize(3);
  #endif
  m_coordMid.resize(4);
  m_coordFace.resize(4);
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

    int compt = 0;
    ENUMERATE_NODE(inode, cell.nodes())
    {
      #if VERIF
      if(m_coord[inode][3] == -123)
      {
        if( m_coord[inode][0] != offset[compt].x() + x ||
            m_coord[inode][1] != offset[compt].y() + y ||
            m_coord[inode][2] != offset[compt].z() + z)
          {
            error() << m_coord[inode][0] << " " << offset[compt].x() + x;
            error() << m_coord[inode][1] << " " << offset[compt].y() + y;
            error() << m_coord[inode][2] << " " << offset[compt].z() + z;
            ARCANE_FATAL("Erreur Node");
          }
      }
      m_coord[inode][3] = -123;
      #endif
      m_coord[inode][0] = offset[compt].x() + x;
      m_coord[inode][1] = offset[compt].y() + y;
      m_coord[inode][2] = offset[compt].z() + z;
      compt++;
      //info() << m_coord[inode][0] << " " << m_coord[inode][1] << " " << m_coord[inode][2];
    }

    ENUMERATE_FACE(iface, cell.faces())
    {
      m_coordMid[iface][0] = offset[compt].x() + x;
      m_coordMid[iface][1] = offset[compt].y() + y;
      m_coordMid[iface][2] = offset[compt].z() + z;
      m_coordMid[iface][3] = offset[compt].b();

      m_coordFace[iface][0] = std::max(std::min(0, faceTupleOffset[compt].x() + x), nx-1); // MinMax pour éviter pos négatives.
      m_coordFace[iface][1] = std::max(std::min(0, faceTupleOffset[compt].y() + y), ny-1);
      m_coordFace[iface][2] = std::max(std::min(0, faceTupleOffset[compt].z() + z), nz-1);

      // Définir les conditions boundary.
      //////////////////// Début conditions boundary /////////////////////////
      qs_vector<MC_Subfacet_Adjacency_Event::Enum> condiBound = getBoundaryCondition(params);
      //Real faceNbr = m_coordFace[iface][0] + nx*(m_coordFace[iface][1] + ny*(m_coordFace[iface][2]));

      // Si la face est au bord du domaine entier.
      if(m_coordFace[iface][0] == x && m_coordFace[iface][1] == y && m_coordFace[iface][2] == z)
      {
        m_boundaryCond[iface] = condiBound[faceTupleOffset[compt].b()];
      }
      else
      {
        Face face = *iface;
        // Si la face est au bord du sous-domaine.
        if(face.isSubDomainBoundary())
        {
          m_boundaryCond[iface] = MC_Subfacet_Adjacency_Event::Transit_Off_Processor;
        }
        
        else
        {
          m_boundaryCond[iface] = MC_Subfacet_Adjacency_Event::Transit_On_Processor;
        }
      }
      //////////////////// Fin conditions boundary /////////////////////////

      compt++;
    }
  }
  //////////////////// Fin Numérotation Node/Face /////////////////////////











  
  int nDomainsPerRank = 1; // SAD set this to 1 for some types of tests
  if( xDom == 0 && yDom == 0 && zDom == 0 )
      if (nRanks == 1)
        nDomainsPerRank = 4;
  


  DecompositionObject ddc(myRank, nRanks, nDomainsPerRank, 0);
  vector<int> myDomainGid = ddc.getAssignedDomainGids(); // Nos sous-domaines
  
  GlobalFccGrid globalGrid(nx, ny, nz, lx, ly, lz);
  
  int nCenters = nRanks*nDomainsPerRank;
  vector<MC_Vector> domainCenter;
  if (xDom == 0 && yDom == 0 && zDom == 0)
      initializeCentersRandomly(nCenters, globalGrid, domainCenter);
  else
      initializeCentersGrid(lx, ly, lz, xDom, yDom, zDom, domainCenter);
  
  qs_assert(domainCenter.size() == nCenters);
  
  // Liste de nos sous-domaines.
  vector<MeshPartition> partition;
  {
      int foremanRank = myRank;
      for (unsigned ii=0; ii<myDomainGid.size(); ++ii)
      {
        partition.push_back(MeshPartition(myDomainGid[ii], ii, foremanRank));
        qs_assert(ddc.getIndex(myDomainGid[ii]) == ii);
      }
  }
  
  CommObject* comm = 0;
  if (nRanks == 1)
      comm = new SharedMemoryCommObject(partition);
  else if (nRanks > 1 && nDomainsPerRank == 1)
      comm = new MpiCommObject(MPI_COMM_WORLD, ddc);
  else
      qs_assert(false);
  
  // S'occupe de lier les cells aux domaines (pas utile ici),
  // et de faire la carte des voisins des sous-domaines (TODO comment faire ça avec arcane ?).
  for (unsigned ii=0; ii<myDomainGid.size(); ++ii)
  {
      if (myRank == 0) { cout << "Building partition " << myDomainGid[ii] << endl; }
      partition[ii].buildMeshPartition(globalGrid, domainCenter, comm);
  }
  
  mpiBarrier(MPI_COMM_WORLD);
  if (myRank ==0 ) { cout << "done building" << endl; }
  mpiBarrier(MPI_COMM_WORLD);
  
  delete comm;
  


  monteCarlo->domain.reserve(myDomainGid.size(),VAR_MEM);
  monteCarlo->domain.Open();
  for (unsigned ii=0; ii<myDomainGid.size(); ++ii)
  {
      if (myRank == 0) { cout << "Building MC_Domain " << ii << endl; }
      monteCarlo->domain.push_back(
        MC_Domain(partition[ii], globalGrid, ddc, params, *monteCarlo->_materialDatabase,
                  params.simulationParams.nGroups));
  }
  monteCarlo->domain.Close();
  
  if (nRanks == 1)
      consistencyCheck(myRank, monteCarlo->domain);
  
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
  monteCarlo->_tallies->InitializeTallies(
      monteCarlo,  
      params.simulationParams.balanceTallyReplications,
      params.simulationParams.fluxTallyReplications,
      params.simulationParams.cellTallyReplications
  );
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

