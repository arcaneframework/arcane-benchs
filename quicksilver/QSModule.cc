// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "QSModule.hh"
#include "arcane/IParallelMng.h"
#include "arcane/IMesh.h"

#include <iostream>
#include "NVTX_Range.hh"
#include "qs_assert.hh"

using namespace Arcane;

bool ordre_qs[]  = {true, true, false, false, false, true};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

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

*/

void QSModule::
initModule()
{
  m_cartesian_mesh = ICartesianMesh::getReference(mesh(), true);
  m_cartesian_mesh->computeDirections();
  initTallies();
  initMesh();
}

void QSModule::
cycleFinalize()
{
  // Update the cumulative tally data.
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
cycleFinalizeTallies()
{
  m_absorb.reduce(Parallel::ReduceSum);
  m_census.reduce(Parallel::ReduceSum);
  m_escape.reduce(Parallel::ReduceSum);
  m_collision.reduce(Parallel::ReduceSum);
  m_end.reduce(Parallel::ReduceSum);
  m_fission.reduce(Parallel::ReduceSum);
  m_produce.reduce(Parallel::ReduceSum);
  m_scatter.reduce(Parallel::ReduceSum);
  m_start.reduce(Parallel::ReduceSum);
  m_source.reduce(Parallel::ReduceSum);
  m_rr.reduce(Parallel::ReduceSum);
  m_split.reduce(Parallel::ReduceSum);
  m_numSegments.reduce(Parallel::ReduceSum);

  info() << "m_start : " << m_start.value();
  info() << "m_source : " << m_source.value();
  info() << "m_rr : " << m_rr.value();
  info() << "m_split : " << m_split.value();
  info() << "m_absorb : " << m_absorb.value();
  info() << "m_scatter : " << m_scatter.value();
  info() << "m_fission : " << m_fission.value();
  info() << "m_produce : " << m_produce.value();
  info() << "m_collision : " << m_collision.value();
  info() << "m_escape : " << m_escape.value();
  info() << "m_census : " << m_census.value();
  info() << "m_numSegments : " << m_numSegments.value();
  info() << "m_end : " << m_end.value();
}

void QSModule::
initMesh()
{
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

  Real lx = options()->getLx();
  Real ly = options()->getLy();
  Real lz = options()->getLz();

  Real dx = lx / m_nx;
  Real dy = ly / m_ny;
  Real dz = lz / m_nz;

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

  UniqueArray<IntegerUniqueArray> faceTupleOffset;
  faceTupleOffset.reserve(6);
  faceTupleOffset.push_back( IntegerUniqueArray{ 0,  0, -1, 5} ); // 13
  faceTupleOffset.push_back( IntegerUniqueArray{-1,  0,  0, 1} ); // 9
  faceTupleOffset.push_back( IntegerUniqueArray{ 0, -1,  0, 3} ); // 11
  faceTupleOffset.push_back( IntegerUniqueArray{ 0,  0,  1, 4} ); // 12
  faceTupleOffset.push_back( IntegerUniqueArray{ 1,  0,  0, 0} ); // 8
  faceTupleOffset.push_back( IntegerUniqueArray{ 0,  1,  0, 2} ); // 10
   

  // Ici, les cells ont déjà une numérotation (cell.uniqueId() (G)  ou  cell.localId() (L))
  // On doit numéroter les nodes selon cell.uniqueId.
  //////////////////// Début Numérotation Node/Face /////////////////////////

  m_coord.resize(4);
  m_coordCm.resize(3);
  m_coordMidCm.resize(3);
  m_coordMid.resize(4);
  m_coordFace.resize(4);
  m_coordCenter.resize(4);

  m_particleCoord.resize(3);
  m_particleVelocity.resize(3);
  m_particleDirCos.resize(3);

  m_total.resize(options()->getNGroups());

  ENUMERATE_CELL(icell, ownCells())
  {
    Cell cell = *icell;
    Int32 unique_id = cell.uniqueId().asInt32();

    Integer index = unique_id;
    Integer x = index % m_nx;
    index /= m_nx;
    Integer y = index % m_ny;
    Integer z = index / m_ny;

    m_coordCenter[icell][0] = 0.0;
    m_coordCenter[icell][1] = 0.0;
    m_coordCenter[icell][2] = 0.0;

    Integer compt = 0;
    ENUMERATE_NODE(inode, cell.nodes())
    {
      m_coord[inode][0] = offset[compt][0] + x;
      m_coord[inode][1] = offset[compt][1] + y;
      m_coord[inode][2] = offset[compt][2] + z;
      m_coord[inode][3] = offset[compt][3];
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

    Real volume = 0;
    MC_Vector cellCenter(m_coordCenter[icell]);

    //info() << "Cell #" << unique_id << " x : " << x << " y : " << y << " z : " << z;

    ENUMERATE_FACE(iface, cell.faces())
    {
      Face face = *iface;

      m_indexArc[iface] = iface.index();

      m_coordMid[iface][0] = offset[compt][0] + x;
      m_coordMid[iface][1] = offset[compt][1] + y;
      m_coordMid[iface][2] = offset[compt][2] + z;
      m_coordMid[iface][3] = offset[compt][3];

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

      Integer compt2 = compt - 8;

      m_coordFace[iface][0] = std::min(std::max(0, faceTupleOffset[compt2][0] + x), (Integer)m_nx-1); // MinMax pour éviter pos négatives.
      m_coordFace[iface][1] = std::min(std::max(0, faceTupleOffset[compt2][1] + y), (Integer)m_ny-1); // TODO : Retirer int cast.
      m_coordFace[iface][2] = std::min(std::max(0, faceTupleOffset[compt2][2] + z), (Integer)m_nz-1);

      //info() << "  Face #" << m_indexArc[iface];

      // Définir les conditions boundary.
      //////////////////// Début conditions boundary /////////////////////////
      UniqueArray<Face_Adjacency_Event> condiBound = getBoundaryCondition();
      //Real faceNbr = m_coordFace[iface][0] + nx*(m_coordFace[iface][1] + ny*(m_coordFace[iface][2]));

      // Si la face est au bord du domaine entier.
      //if(m_coordFace[iface][0] == x && m_coordFace[iface][1] == y && m_coordFace[iface][2] == z)
      if(face.isSubDomainBoundary())
      {
        m_boundaryCond[iface] = condiBound[faceTupleOffset[compt2][3]];
        //info() << "    Bord du domaine entier " << face.nbCell();

      }
      // Si la face est au bord du sous-domaine.
      else if(face.frontCell().owner() != face.backCell().owner())
      {
        // info() << "    Bord du sous-domaine " << face.nbCell();
        // info() <<  "x : " << x << " y: " << y << " z: " << z 
        // << " xx : " << faceTupleOffset[compt2].x() 
        // << " yy : " << faceTupleOffset[compt2].y() 
        // << " zz : " << faceTupleOffset[compt2].z() 
        // << " compt : " << compt2;
        m_boundaryCond[iface] = Face_Adjacency_Event::Transit_Off_Processor;
      }
      
      else
      {
        //info() << "    Intra " << face.nbCell();
        m_boundaryCond[iface] = Face_Adjacency_Event::Transit_On_Processor;
      }
      
      //////////////////// Fin conditions boundary /////////////////////////

      //////////////////// Début Volume Cell /////////////////////////


      for (Integer i = 0; i < 4; i++)
      {
        Integer first_pos_node = (ordre_qs[iface.index()] ? ((i == 3) ? 0 : i+1) : i);
        Integer second_pos_node = (ordre_qs[iface.index()] ? i : ((i == 3) ? 0 : i+1));

        Node first_node = face.node(first_pos_node);
        Node second_node = face.node(second_pos_node);

        MC_Vector aa = MC_Vector(m_coordCm[first_node]) - cellCenter;
        MC_Vector bb = MC_Vector(m_coordCm[second_node]) - cellCenter;
        MC_Vector cc = MC_Vector(m_coordMidCm[iface]) - cellCenter;

        volume += abs(aa.Dot(bb.Cross(cc)));
      }

      
      compt++;
    }
    volume /= 6.0;

    m_volume[cell] = volume;
    //////////////////// Fin Volume Cell /////////////////////////

    // for(Integer i = 0; i < options()->getNGroups(); i++)
    // {
    //   m_total[icell][i] = 0.0;
    // }

    //m_cellNumberDensity[icell] = 1.0;
    //m_sourceTally[icell] = 0;

    
  }
  //////////////////// Fin Numérotation Node/Face /////////////////////////

  m_total.fill(0.0);
  m_cellNumberDensity.fill(1.0);
  m_sourceTally.fill(0);

  info() << "Finished initMesh";
}

void QSModule::
initTallies()
{
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

  Integer sizeOfSFT = options()->getNGroups();
  m_scalarFluxTally.resize(sizeOfSFT);

  ENUMERATE_CELL(icell, ownCells())
  {
    m_cellTally[icell] = 0.0;
    for (Integer i = 0; i < sizeOfSFT; i++)
    {
      m_scalarFluxTally[icell][i] = 0.0;
      // TODO : Voir pour ajouter atomic.
    }
  }
}

UniqueArray<Face_Adjacency_Event> QSModule::
getBoundaryCondition()
{
  UniqueArray<Face_Adjacency_Event> bc(6);
  if (options()->getBoundaryCondition() == eBoundaryCondition::reflect)
  {
    bc = UniqueArray<Face_Adjacency_Event>(6, Face_Adjacency_Event::Boundary_Reflection);
  }
  else if (options()->getBoundaryCondition() == eBoundaryCondition::escape)
  {
    bc = UniqueArray<Face_Adjacency_Event>(6, Face_Adjacency_Event::Boundary_Escape);
  }
  else if  (options()->getBoundaryCondition() == eBoundaryCondition::octant)
  {
    for (unsigned ii=0; ii<6; ++ii)
    {
      if (ii % 2 == 0) bc[ii] = Face_Adjacency_Event::Boundary_Escape;
      if (ii % 2 == 1) bc[ii] = Face_Adjacency_Event::Boundary_Reflection;
    }
  }
  else
  {
    qs_assert(false);
  }
  return bc;
}
