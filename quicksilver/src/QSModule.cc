// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* QSModule.cc                                                 (C) 2000-2022 */
/*                                                                           */
/* Module principal QAMA                                                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "QSModule.hh"
#include <iostream>

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

  VariableNodeReal3& node_coord = nodesCoordinates();

  m_global_deltat = options()->getDt();

  m_e_min = options()->getEMin();
  m_e_max = options()->getEMax();
  m_n_groups = options()->getNGroups();

  m_lx = options()->getLx(); // TODO : Récupérer directement les <length> du cartesian mesh.
  m_ly = options()->getLy();
  m_lz = options()->getLz();

  // Voir si la parallélisation ici sert à quelque chose
  // sachant qu'il faut gérer la répartition des faces par threads.
  ENUMERATE_CELL (icell, ownCells()) {
    Cell cell = *icell;

    m_cell_center_coord[icell][MD_DirX] = 0.0;
    m_cell_center_coord[icell][MD_DirY] = 0.0;
    m_cell_center_coord[icell][MD_DirZ] = 0.0;

    ENUMERATE_FACE (iface, cell.faces()) {
      Face face = *iface;

      m_face_center_coord[iface][MD_DirX] = 0.0;
      m_face_center_coord[iface][MD_DirY] = 0.0;
      m_face_center_coord[iface][MD_DirZ] = 0.0;

      ENUMERATE_NODE (inode, face.nodes()) {
        m_face_center_coord[iface][MD_DirX] += node_coord[inode][MD_DirX];
        m_face_center_coord[iface][MD_DirY] += node_coord[inode][MD_DirY];
        m_face_center_coord[iface][MD_DirZ] += node_coord[inode][MD_DirZ];
      }

      m_face_center_coord[iface][MD_DirX] /= 4;
      m_face_center_coord[iface][MD_DirY] /= 4;
      m_face_center_coord[iface][MD_DirZ] /= 4;
      
      m_cell_center_coord[icell][MD_DirX] += m_face_center_coord[iface][MD_DirX];
      m_cell_center_coord[icell][MD_DirY] += m_face_center_coord[iface][MD_DirY];
      m_cell_center_coord[icell][MD_DirZ] += m_face_center_coord[iface][MD_DirZ];
    }

    m_cell_center_coord[icell][MD_DirX] /= 6;
    m_cell_center_coord[icell][MD_DirY] /= 6;
    m_cell_center_coord[icell][MD_DirZ] /= 6;

    Real volume = 0;
    ENUMERATE_FACE (iface, cell.faces()) {
      Face face = *iface;

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

        Real3 aa = node_coord[first_node] - m_cell_center_coord[icell];
        Real3 bb = node_coord[second_node] - m_cell_center_coord[icell];
        Real3 cc = m_face_center_coord[iface] - m_cell_center_coord[icell];

        volume += math::abs(math::dot(aa, math::cross(bb, cc)));
      }
    }
    volume /= 6.0;

    m_volume[icell] = volume;
  }

  m_normalFace.resize(6);
  ENUMERATE_FACE(iface, ((CellVectorView) ownCells().view())[0].faces())
  {
    Face face = (*iface);
    Real3 n0 = m_face_center_coord[iface];
    Real3 n1 = node_coord[face.node(0)];
    Real3 n2 = node_coord[face.node(1)];

    // Les trois premières faces sont les faces externes,
    // les trois suivantes sont internes 
    // (en regardant la cellule (0, 0, 0)).
    Real sens = (iface.index() < 3 ? -1.0 : 1.0);

    // Calcul du vecteur normal à la face, puis normalisation.
    Real aa = (((n1.y - n0.y) * (n2.z - n0.z)) - ((n1.z - n0.z) * (n2.y - n0.y)) == 0 ? 0 : sens);
    Real bb = (((n1.z - n0.z) * (n2.x - n0.x)) - ((n1.x - n0.x) * (n2.z - n0.z)) == 0 ? 0 : sens);
    Real cc = (((n1.x - n0.x) * (n2.y - n0.y)) - ((n1.y - n0.y) * (n2.x - n0.x)) == 0 ? 0 : sens);

    m_normalFace[iface.index()] = Real3(aa, bb, cc);
  }

  m_total.resize(m_n_groups());
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
getBoundaryCondition(const Integer& pos)
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