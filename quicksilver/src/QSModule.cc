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
 * au maillage.
 */
void QSModule::
initModule()
{
  initMesh();

  // Initialisation de la sortie CSV.
  ISimpleOutput* csv = ServiceBuilder<ISimpleOutput>(subDomain()).getSingleton();
  csv->init("QAMA", ";");

  // On crée les lignes "Proc" ici pour les avoir en haut du csv du proc0.
  csv->addRow("Sampling duration (Proc)");
  csv->addRow("Tracking duration (Proc)");
  csv->addRow("m_start (Proc)");
  csv->addRow("m_source (Proc)");
  csv->addRow("m_rr (Proc)");
  csv->addRow("m_split (Proc)");
  csv->addRow("m_absorb (Proc)");
  csv->addRow("m_scatter (Proc)");
  csv->addRow("m_fission (Proc)");
  csv->addRow("m_produce (Proc)");
  csv->addRow("m_collision (Proc)");
  csv->addRow("m_escape (Proc)");
  csv->addRow("m_census (Proc)");
  csv->addRow("m_num_segments (Proc)");
  csv->addRow("m_end (Proc)");
  csv->addRow("m_incoming (Proc)");
  csv->addRow("m_outgoing (Proc)");
  csv->addRow("sum_scalar_flux_tally (Proc)");
}

/**
 * @brief Méthode permettant d'afficher les informations de fin d'itération.
 */
void QSModule::
cycleFinalize()
{
  info() << "End iteration #" << m_global_iteration();
  if(parallelMng()->commSize() == 1) info() << "  Informations:";
  else  info() << "  Informations:                                            (variable name): Sum, [Min, Max, Avg]";

  if (m_global_iteration() == options()->getNSteps())
    subDomain()->timeLoopMng()->stopComputeLoop(true);
}

/**
 * @brief Méthode appelée à la fin de la boucle en temps.
 */
void QSModule::
endModule()
{
  if(options()->getCsvFile() != "") {
    ISimpleOutput* csv = ServiceBuilder<ISimpleOutput>(subDomain()).getSingleton();
    //csv->print();
    csv->writeFile(options()->getCsvFile());
  }
}

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

  ICartesianMesh* cartesian_mesh = ICartesianMesh::getReference(mesh(), true);
  cartesian_mesh->computeDirections();

  // On recherche le nombre de cellules en x, y, z.
  Int64 m_nx, m_ny, m_nz;
  {
    CellDirectionMng cdm(cartesian_mesh->cellDirection(MD_DirX));
    m_nx = cdm.globalNbCell();
  }
  {
    CellDirectionMng cdm(cartesian_mesh->cellDirection(MD_DirY));
    m_ny = cdm.globalNbCell();
  }
  {
    CellDirectionMng cdm(cartesian_mesh->cellDirection(MD_DirZ));
    m_nz = cdm.globalNbCell();
  }

  VariableNodeReal3& node_coord = nodesCoordinates();

  ICartesianMeshGenerationInfo* aaa = ICartesianMeshGenerationInfo::getReference(mesh(), false);

  Real3 lenght = aaa->globalLength();
  m_lx = lenght[MD_DirX];
  m_ly = lenght[MD_DirY];
  m_lz = lenght[MD_DirZ];

  m_global_deltat = options()->getDt();

  m_e_min = options()->getEMin();
  m_e_max = options()->getEMax();
  m_n_groups = options()->getNGroups();

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

  m_normal_face.resize(6);
  ENUMERATE_FACE (iface, ((CellVectorView)ownCells().view())[0].faces()) {
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

    m_normal_face[iface.index()] = Real3(aa, bb, cc);
  }

  m_total.resize(m_n_groups());
  m_total.fill(0.0);
  m_cell_number_density.fill(1.0);
  m_source_tally.fill(0);
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