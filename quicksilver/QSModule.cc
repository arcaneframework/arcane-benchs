// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "QSModule.hh"
#include <iostream>
#include <arcane/Concurrency.h>
#include "NVTX_Range.hh"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Méthode permettant de lancer l'initialisation des grandeurs
 * au maillage et des tallies.
 */
void QSModule::
initModule()
{
  m_cartesian_mesh = ICartesianMesh::getReference(mesh(), true);
  m_cartesian_mesh->computeDirections();
  initMesh();
  initTallies();
}

/**
 * @brief Méthode permettant d'afficher les informations de fin d'itération.
 */
void QSModule::
cycleFinalize()
{
  cycleFinalizeTallies();

  if (m_global_iteration() == options()->getNSteps())
    subDomain()->timeLoopMng()->stopComputeLoop(true);
}

/**
 * @brief Méthode appelée à la fin de la boucle en temps.
 */
void QSModule::
endModule() {}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Méthode permettant d'initialiser les grandeurs au maillage (entre
 * autres choses).
 */
void QSModule::
initMesh()
{
  info() << "Initialisation des grandeurs/variables";

  // On recherche le nombre de cellules en x, y, z.
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

  m_e_min = options()->getEMin();
  m_e_max = options()->getEMax();
  m_n_groups = options()->getNGroups();

  m_lx = options()->getLx();
  m_ly = options()->getLy();
  m_lz = options()->getLz();

  Real dx = m_lx() / m_nx;
  Real dy = m_ly() / m_ny;
  Real dz = m_lz() / m_nz;

  // Extrait de GlobalFccGrid.cc.
  // Permet de retrouver la position de chaque noeud dans l'espace.
  UniqueArray<IntegerUniqueArray> offset;
  offset.reserve(14);
  offset.push_back(IntegerUniqueArray{ 0, 0, 0, 0 }); // 0
  offset.push_back(IntegerUniqueArray{ 1, 0, 0, 0 }); // 1
  offset.push_back(IntegerUniqueArray{ 1, 1, 0, 0 }); // 3
  offset.push_back(IntegerUniqueArray{ 0, 1, 0, 0 }); // 2

  offset.push_back(IntegerUniqueArray{ 0, 0, 1, 0 }); // 4
  offset.push_back(IntegerUniqueArray{ 1, 0, 1, 0 }); // 5
  offset.push_back(IntegerUniqueArray{ 1, 1, 1, 0 }); // 7
  offset.push_back(IntegerUniqueArray{ 0, 1, 1, 0 }); // 6

  offset.push_back(IntegerUniqueArray{ 0, 0, 0, 3 }); // 13
  offset.push_back(IntegerUniqueArray{ 0, 0, 0, 1 }); // 9
  offset.push_back(IntegerUniqueArray{ 0, 0, 0, 2 }); // 11
  offset.push_back(IntegerUniqueArray{ 0, 0, 1, 3 }); // 12
  offset.push_back(IntegerUniqueArray{ 1, 0, 0, 1 }); // 8
  offset.push_back(IntegerUniqueArray{ 0, 1, 0, 2 }); // 10

  m_coord_cm.resize(3);
  m_coord_mid_cm.resize(3);
  m_coord_center.resize(4);
  m_total.resize(m_n_groups());

  Arcane::ParallelLoopOptions options;

  // Exécute la boucle par parties d'environ 50 mailles.
  // options.setGrainSize(50);
  // options.setMaxThread(3);

  arcaneParallelForeach(ownCells(), options, [&](CellVectorView cells) {
    ENUMERATE_CELL (icell, cells) {
      Cell cell = *icell;

      // Position de la cellule dans le maillage entier.
      Integer index = cell.uniqueId().asInt32();
      Integer x = index % m_nx;
      index /= m_nx;
      Integer y = index % m_ny;
      Integer z = index / m_ny;

      m_coord_center[icell][MD_DirX] = 0.0;
      m_coord_center[icell][MD_DirY] = 0.0;
      m_coord_center[icell][MD_DirZ] = 0.0;

      ENUMERATE_NODE (inode, cell.nodes()) {
        Real m_coordX = offset[inode.index()][0] + x;
        Real m_coordY = offset[inode.index()][1] + y;
        Real m_coordZ = offset[inode.index()][2] + z;

        // Coordonnées du noeud en cm.
        m_coord_cm[inode][MD_DirX] = m_coordX * dx;
        m_coord_cm[inode][MD_DirY] = m_coordY * dy;
        m_coord_cm[inode][MD_DirZ] = m_coordZ * dz;

        m_coord_center[icell][MD_DirX] += m_coord_cm[inode][MD_DirX];
        m_coord_center[icell][MD_DirY] += m_coord_cm[inode][MD_DirY];
        m_coord_center[icell][MD_DirZ] += m_coord_cm[inode][MD_DirZ];
      }

      m_coord_center[icell][MD_DirX] /= cell.nbNode();
      m_coord_center[icell][MD_DirY] /= cell.nbNode();
      m_coord_center[icell][MD_DirZ] /= cell.nbNode();

      Integer compt = 8;
      Real volume = 0;
      Real3 cellCenter = avToReal3(m_coord_center[icell]);

      ENUMERATE_FACE (iface, cell.faces()) {
        Face face = *iface;

        m_index_arc[iface] = iface.index();

        Real m_coordX = offset[compt][0] + x;
        Real m_coordY = offset[compt][1] + y;
        Real m_coordZ = offset[compt][2] + z;
        Real m_coordB = offset[compt][3];

        // Coordonnées du millieux de la face (point commun à toutes les facets
        // de la face).
        m_coord_mid_cm[iface][MD_DirX] = m_coordX * dx;
        m_coord_mid_cm[iface][MD_DirY] = m_coordY * dy;
        m_coord_mid_cm[iface][MD_DirZ] = m_coordZ * dz;

        if (m_coordB == 1) {
          m_coord_mid_cm[iface][MD_DirY] += dy / 2;
          m_coord_mid_cm[iface][MD_DirZ] += dz / 2;
        }
        else if (m_coordB == 2) {
          m_coord_mid_cm[iface][MD_DirX] += dx / 2;
          m_coord_mid_cm[iface][MD_DirZ] += dz / 2;
        }
        else {
          m_coord_mid_cm[iface][MD_DirX] += dx / 2;
          m_coord_mid_cm[iface][MD_DirY] += dy / 2;
        }

        // info() << "  Face #" << m_index_arc[iface];

        // Si la face est au bord du domaine entier.
        if (face.isSubDomainBoundary()) {
          // D'origine, dans le code QS, dans le cas octant :
          //   les faces 0, 1, 2 sont escape 
          //   les faces 3, 4, 5 sont reflection
          m_boundary_cond[iface] = getBoundaryCondition(iface.index());
        }
        // Face interne au sous-domaine ou au bord du sous-domaine.
        else {
          m_boundary_cond[iface] = ParticleEvent::cellChange;
        }

        for (Integer i = 0; i < 4; i++) {
          Node first_node = face.node(i);
          Node second_node = face.node(((i == 3) ? 0 : i + 1));

          Real3 aa = avToReal3(m_coord_cm[first_node]) - cellCenter;
          Real3 bb = avToReal3(m_coord_cm[second_node]) - cellCenter;
          Real3 cc = avToReal3(m_coord_mid_cm[iface]) - cellCenter;

          volume += math::abs(math::dot(aa, math::cross(bb, cc)));
        }
        compt++;
      }
      volume /= 6.0;

      m_volume[cell] = volume;
    }
  });

  m_total.fill(0.0);
  m_cell_number_density.fill(1.0);
  m_source_tally.fill(0);
}

/**
 * @brief Méthode permettant d'initialiser les tallies.
 */
void QSModule::
initTallies()
{
  info() << "Init tallies";
  m_absorb = 0; // Number of particles absorbed
  m_census = 0; // Number of particles that enter census
  m_escape = 0; // Number of particles that escape
  m_collision = 0; // Number of collosions
  m_end = 0; // Number of particles at end of cycle
  m_fission = 0; // Number of fission events
  m_produce = 0; // Number of particles created by collisions
  m_scatter = 0; // Number of scatters
  m_start = 0; // Number of particles at beginning of cycle
  m_source = 0; // Number of particles sourced in
  m_rr = 0; // Number of particles Russian Rouletted in population control
  m_split = 0; // Number of particles split in population control
  m_num_segments = 0; // Number of segements

  m_scalar_flux_tally.resize(m_n_groups());
  m_scalar_flux_tally.fill(0.0);
}

/**
 * @brief Méthode permettant de récupérer les tallies de tous les sous-domaines
 * et de les afficher.
 */
void QSModule::
cycleFinalizeTallies()
{

  // Somme des m_scalar_flux_tally.
  Real sum_scalar_flux_tally = 0.0;
  ENUMERATE_CELL(icell, ownCells()){
    for(Integer i = 0; i < m_n_groups(); i++){
      sum_scalar_flux_tally += m_scalar_flux_tally[icell][i];
      m_scalar_flux_tally[icell][i] = 0.0;
    }
  }

  m_absorb.reduce(Parallel::ReduceSum);
  m_census.reduce(Parallel::ReduceSum);
  m_escape.reduce(Parallel::ReduceSum);
  m_collision.reduce(Parallel::ReduceSum);
  m_fission.reduce(Parallel::ReduceSum);
  m_produce.reduce(Parallel::ReduceSum);
  m_scatter.reduce(Parallel::ReduceSum);
  m_num_segments.reduce(Parallel::ReduceSum);
  m_source.reduce(Parallel::ReduceSum);
  m_rr.reduce(Parallel::ReduceSum);
  m_split.reduce(Parallel::ReduceSum);
  m_start.reduce(Parallel::ReduceSum);
  m_end.reduce(Parallel::ReduceSum);
  sum_scalar_flux_tally = mesh()->parallelMng()->reduce(Parallel::ReduceSum, sum_scalar_flux_tally);

  info() << "End iteration #" << m_global_iteration();
  info() << "  Informations:";
  info() << "    Number of particles at beginning of cycle                    "
            "(m_start): "
         << m_start.value();
  info() << "    Number of particles sourced in                              "
            "(m_source): "
         << m_source.value();
  info() << "    Number of particles Russian Rouletted in population control   "
            "  (m_rr): "
         << m_rr.value();
  info() << "    Number of particles split in population control              "
            "(m_split): "
         << m_split.value();
  info() << "    Number of particles absorbed                                "
            "(m_absorb): "
         << m_absorb.value();
  info() << "    Number of scatters                                         "
            "(m_scatter): "
         << m_scatter.value();
  info() << "    Number of fission events                                   "
            "(m_fission): "
         << m_fission.value();
  info() << "    Number of particles created by collisions                  "
            "(m_produce): "
         << m_produce.value();
  info() << "    Number of collisions                                     "
            "(m_collision): "
         << m_collision.value();
  info() << "    Number of particles that escape                             "
            "(m_escape): "
         << m_escape.value();
  info() << "    Number of particles that enter census                       "
            "(m_census): "
         << m_census.value();
  info() << "    Number of segements                                   "
            "(m_num_segments): "
         << m_num_segments.value();
  info() << "    Number of particles at end of cycle                           "
            " (m_end): "
         << m_end.value();
  info() << "    Particles contribution to the scalar flux     "
            " (sum_scalar_flux_tally): "
         << sum_scalar_flux_tally;

  m_start = 0;
  m_source = 0;
  m_rr = 0;
  m_split = 0;

  m_absorb = 0;
  m_census = 0;
  m_escape = 0;
  m_collision = 0;
  m_fission = 0;
  m_produce = 0;
  m_scatter = 0;
  m_num_segments = 0;
  m_end = 0;
}

/**
 * @brief Méthode permettant de récupérer la condition aux bords de maillage.
 *
 * @param pos Si la condition est 'octant', on choisi soit escape si pos pair,
 * soit reflection sinon.
 * @return ParticleEvent La condition aux bords de maillage.
 */
ParticleEvent QSModule::
getBoundaryCondition(Integer pos)
{
  switch (options()->getBoundaryCondition()) {
  case eBoundaryCondition::REFLECT:
    return ParticleEvent::reflection;

  case eBoundaryCondition::ESCAPE:
    return ParticleEvent::escape;

  case eBoundaryCondition::OCTANT:
    if (pos < 3)
      return ParticleEvent::escape;
    else
      return ParticleEvent::reflection;

  default:
    ARCANE_ASSERT(false, "Boundary condition undefined.");
    return ParticleEvent::undefined;
  }
}