// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "QSModule.hh"
#include <iostream>
#include "NVTX_Range.hh"
#include "qs_assert.hh"

using namespace Arcane;


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void QSModule::
initModule()
{
  m_cartesian_mesh = ICartesianMesh::getReference(mesh(), true);
  m_cartesian_mesh->computeDirections();
  initMesh();
  initTallies();
}

void QSModule::
cycleFinalize()
{
  cycleFinalizeTallies();

  if(m_global_iteration() == options()->getNSteps())
    subDomain()->timeLoopMng()->stopComputeLoop(true);
}

void QSModule::
endModule()
{

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void QSModule::
initMesh()
{
  info() << "Initialisation des grandeurs/variables";

  Int64 m_nx, m_ny, m_nz;
  {
    CellDirectionMng cdm(m_cartesian_mesh->cellDirection(MD_DirX));
    m_nx = cdm.globalNbCell();
  }
  {
    CellDirectionMng cdm(m_cartesian_mesh->cellDirection(MD_DirY));
    m_ny = cdm.globalNbCell();
  }
  {
    CellDirectionMng cdm(m_cartesian_mesh->cellDirection(MD_DirZ));
    m_nz = cdm.globalNbCell();
  }

  m_global_deltat = options()->getDt();

  m_eMin = options()->getEMin();
  m_eMax = options()->getEMax();
  m_nGroups = options()->getNGroups();

  m_lx = options()->getLx();
  m_ly = options()->getLy();
  m_lz = options()->getLz();

  Real dx = m_lx() / m_nx;
  Real dy = m_ly() / m_ny;
  Real dz = m_lz() / m_nz;

  // Extrait de GlobalFccGrid.cc.
  UniqueArray<IntegerUniqueArray> offset;
  offset.reserve(14);
  offset.push_back(IntegerUniqueArray{0, 0, 0, 0}); // 0
  offset.push_back(IntegerUniqueArray{1, 0, 0, 0}); // 1
  offset.push_back(IntegerUniqueArray{1, 1, 0, 0}); // 3
  offset.push_back(IntegerUniqueArray{0, 1, 0, 0}); // 2

  offset.push_back(IntegerUniqueArray{0, 0, 1, 0}); // 4
  offset.push_back(IntegerUniqueArray{1, 0, 1, 0}); // 5
  offset.push_back(IntegerUniqueArray{1, 1, 1, 0}); // 7
  offset.push_back(IntegerUniqueArray{0, 1, 1, 0}); // 6

  offset.push_back(IntegerUniqueArray{0, 0, 0, 3}); // 13
  offset.push_back(IntegerUniqueArray{0, 0, 0, 1}); // 9
  offset.push_back(IntegerUniqueArray{0, 0, 0, 2}); // 11
  offset.push_back(IntegerUniqueArray{0, 0, 1, 3}); // 12
  offset.push_back(IntegerUniqueArray{1, 0, 0, 1}); // 8
  offset.push_back(IntegerUniqueArray{0, 1, 0, 2}); // 10
  
  // Ici, les cells ont déjà une numérotation (cell.uniqueId() (G)  ou  cell.localId() (L))
  // On doit numéroter les nodes selon cell.uniqueId.
  //////////////////// Début Numérotation Node/Face /////////////////////////

  m_coordCm.resize(3);
  m_coordMidCm.resize(3);
  m_coordCenter.resize(4);

  m_total.resize(options()->getNGroups());

  ENUMERATE_CELL(icell, ownCells())
  {
    Cell cell = *icell;

    Integer index = cell.uniqueId().asInt32();
    Integer x = index % m_nx;
    index /= m_nx;
    Integer y = index % m_ny;
    Integer z = index / m_ny;

    m_coordCenter[icell][MD_DirX] = 0.0;
    m_coordCenter[icell][MD_DirY] = 0.0;
    m_coordCenter[icell][MD_DirZ] = 0.0;

    ENUMERATE_NODE(inode, cell.nodes())
    {
      Real m_coordX = offset[inode.index()][0] + x;
      Real m_coordY = offset[inode.index()][1] + y;
      Real m_coordZ = offset[inode.index()][2] + z;

      m_coordCm[inode][MD_DirX] = m_coordX * dx;
      m_coordCm[inode][MD_DirY] = m_coordY * dy;
      m_coordCm[inode][MD_DirZ] = m_coordZ * dz;

      m_coordCenter[icell][MD_DirX] += m_coordCm[inode][MD_DirX];
      m_coordCenter[icell][MD_DirY] += m_coordCm[inode][MD_DirY];
      m_coordCenter[icell][MD_DirZ] += m_coordCm[inode][MD_DirZ];
    }
    m_coordCenter[icell][MD_DirX] /= cell.nbNode();
    m_coordCenter[icell][MD_DirY] /= cell.nbNode();
    m_coordCenter[icell][MD_DirZ] /= cell.nbNode();

    Integer compt = 8;
    Real volume = 0;
    MC_Vector cellCenter(m_coordCenter[icell]);

    ENUMERATE_FACE(iface, cell.faces())
    {
      Face face = *iface;

      m_indexArc[iface] = iface.index();

      Real m_coordX = offset[compt][0] + x;
      Real m_coordY = offset[compt][1] + y;
      Real m_coordZ = offset[compt][2] + z;
      Real m_coordB = offset[compt][3];

      m_coordMidCm[iface][MD_DirX] = m_coordX * dx;
      m_coordMidCm[iface][MD_DirY] = m_coordY * dy;
      m_coordMidCm[iface][MD_DirZ] = m_coordZ * dz;

      if(m_coordB == 1)
      {
        m_coordMidCm[iface][MD_DirY] += dy / 2;
        m_coordMidCm[iface][MD_DirZ] += dz / 2;
      }
      else if(m_coordB == 2)
      {
        m_coordMidCm[iface][MD_DirX] += dx / 2;
        m_coordMidCm[iface][MD_DirZ] += dz / 2;
      }
      else
      {
        m_coordMidCm[iface][MD_DirX] += dx / 2;
        m_coordMidCm[iface][MD_DirY] += dy / 2;
      }

      //info() << "  Face #" << m_indexArc[iface];

      // Si la face est au bord du domaine entier.
      if(face.isSubDomainBoundary())
      {
        // D'origine, dans le cas octant : 
        //   les faces 0, 1, 2 sont escape (donc compt = 8, 9, 10 => getBoundaryCondition(impair))
        //   les faces 3, 4, 5 sont reflection (donc compt = 11, 12, 13 => getBoundaryCondition(pair))
        m_boundaryCond[iface] = getBoundaryCondition((compt/11) + 1); 
      }
      // Si la face est au bord du sous-domaine.
      else if(face.frontCell().owner() != face.backCell().owner())
      {
        m_boundaryCond[iface] = particleEvent::subDChange;
      }
      // Face interne au sous-domaine.
      else
      {
        m_boundaryCond[iface] = particleEvent::cellChange;
      }
      
      for (Integer i = 0; i < 4; i++)
      {
        Node first_node  = face.node(i);
        Node second_node = face.node(((i == 3) ? 0 : i+1));

        MC_Vector aa = MC_Vector(m_coordCm[first_node]) - cellCenter;
        MC_Vector bb = MC_Vector(m_coordCm[second_node]) - cellCenter;
        MC_Vector cc = MC_Vector(m_coordMidCm[iface]) - cellCenter;

        volume += abs(aa.Dot(bb.Cross(cc)));
      }
      compt++;
    }
    volume /= 6.0;

    m_volume[cell] = volume;
  }

  m_total.fill(0.0);
  m_cellNumberDensity.fill(1.0);
  m_sourceTally.fill(0);
}

void QSModule::
initTallies()
{
  info() << "Init tallies";
  m_absorb = 0;      // Number of particles absorbed
  m_census = 0;      // Number of particles that enter census
  m_escape = 0;      // Number of particles that escape
  m_collision = 0;   // Number of collosions
  m_end = 0;         // Number of particles at end of cycle
  m_fission = 0;     // Number of fission events
  m_produce = 0;     // Number of particles created by collisions
  m_scatter = 0;     // Number of scatters
  m_start = 0;       // Number of particles at beginning of cycle
  m_source = 0;      // Number of particles sourced in
  m_rr = 0;          // Number of particles Russian Rouletted in population control
  m_split = 0;       // Number of particles split in population control
  m_numSegments = 0; // Number of segements

  m_scalarFluxTally.resize(options()->getNGroups());
  m_scalarFluxTally.fill(0.0);
  m_cellTally.fill(0.0);
}

void QSModule::
cycleFinalizeTallies()
{
  m_absorb.reduce(Parallel::ReduceSum);
  m_census.reduce(Parallel::ReduceSum);
  m_escape.reduce(Parallel::ReduceSum);
  m_collision.reduce(Parallel::ReduceSum);
  m_fission.reduce(Parallel::ReduceSum);
  m_produce.reduce(Parallel::ReduceSum);
  m_scatter.reduce(Parallel::ReduceSum);
  m_numSegments.reduce(Parallel::ReduceSum);
  m_source.reduce(Parallel::ReduceSum);
  m_rr.reduce(Parallel::ReduceSum);
  m_split.reduce(Parallel::ReduceSum);
  m_start.reduce(Parallel::ReduceSum); // Non utilisée dans code d'origine.
  m_end.reduce(Parallel::ReduceSum); // Non utilisée dans code d'origine.

  info() << "End iteration #" << m_global_iteration();
  info() << "  Informations:";
  //info() << "Number of particles at beginning of cycle (m_start): " << m_start.value();
  info() << "    Number of particles sourced in                          (m_source): " << m_source.value();
  info() << "    Number of particles Russian Rouletted in population control (m_rr): " << m_rr.value();
  info() << "    Number of particles split in population control          (m_split): " << m_split.value();
  info() << "    Number of particles absorbed                            (m_absorb): " << m_absorb.value();
  info() << "    Number of scatters                                     (m_scatter): " << m_scatter.value();
  info() << "    Number of fission events                               (m_fission): " << m_fission.value();
  info() << "    Number of particles created by collisions              (m_produce): " << m_produce.value();
  info() << "    Number of collisions                                 (m_collision): " << m_collision.value();
  info() << "    Number of particles that escape                         (m_escape): " << m_escape.value();
  info() << "    Number of particles that enter census                   (m_census): " << m_census.value();
  info() << "    Number of segements                                (m_numSegments): " << m_numSegments.value();
  //info() << "Number of particles at end of cycle (m_end): " << m_end.value();

  m_absorb = 0;
  m_census = 0;
  m_escape = 0;
  m_collision = 0;
  m_fission = 0;
  m_produce = 0;
  m_scatter = 0;
  m_numSegments = 0;

  m_source = 0;
  m_rr = 0;
  m_split = 0;
}


particleEvent QSModule::
getBoundaryCondition(Integer pos)
{
  switch (options()->getBoundaryCondition())
  {
  case eBoundaryCondition::REFLECT:
    return particleEvent::reflection;

  case eBoundaryCondition::ESCAPE:
    return particleEvent::escape;

  case eBoundaryCondition::OCTANT:
    if (pos % 2 == 0) return particleEvent::escape;
    if (pos % 2 == 1) return particleEvent::reflection;
  
  default:
    qs_assert(false);
    return particleEvent::undefined;
  }
}


/*
1
*I-QS         particle_count : 1017740
*I-QS         Module Quicksilver cycleFinalize
*I-QS         m_start : 0
*I-QS         m_source : 99968
*I-QS         m_rr : 0
*I-QS         m_split : 900314
*I-QS         m_absorb : 35322
*I-QS         m_scatter : 882215
*I-QS         m_fission : 44173
*I-QS         m_produce : 52902
*I-QS         m_collision : 961710
*I-QS         m_escape : 0
*I-QS         m_census : 973689
*I-QS         m_numSegments : 1955279
*I-QS         m_end : 0

2
*I-QS         m_start : 0
*I-QS         m_source : 99968
*I-QS         m_rr : 0
*I-QS         m_split : 0
*I-QS         m_absorb : 173905
*I-QS         m_scatter : 4336085
*I-QS         m_fission : 216560
*I-QS         m_produce : 260030
*I-QS         m_collision : 4726550
*I-QS         m_escape : 0
*I-QS         m_census : 943222
*I-QS         m_numSegments : 5835735
*I-QS         m_end : 0

*I-TrackingMC P0 - End SubIter #6 : Total number of particles processed : 1305004
*I-QS         End iteration #1
*I-QS           Informations:
*I-QS             Number of particles sourced in                          (m_source): 100000
*I-QS             Number of particles Russian Rouletted in population control (m_rr): 0
*I-QS             Number of particles split in population control          (m_split): 900000
*I-QS             Number of particles absorbed                            (m_absorb): 203537
*I-QS             Number of scatters                                     (m_scatter): 5084900
*I-QS             Number of fission events                               (m_fission): 253874
*I-QS             Number of particles created by collisions              (m_produce): 406376
*I-QS             Number of collisions                                 (m_collision): 5542311
*I-QS             Number of particles that escape                         (m_escape): 0
*I-QS             Number of particles that enter census                   (m_census): 948965
*I-QS             Number of segements                                (m_numSegments): 7042210




*/