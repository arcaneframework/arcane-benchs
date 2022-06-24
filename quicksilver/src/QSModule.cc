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
  mesh()->modifier()->setDynamic(true);

  if(options()->getLoadBalancingMat() || options()->getLoadBalancingLoop() > 0) {
    m_criterion_lb_cell.fill(1.);
    m_criterion_lb_face.fill(1.);
  }

  initMesh();
  m_csv = ServiceBuilder<ISimpleTableOutput>(subDomain()).getSingleton();

  // Initialisation de la sortie CSV.
  if(options()->getCsvName() != "")
    m_csv->init(options()->getCsvName(), ";");
  else
    m_csv->init("QAMA_P@proc_id@", ";");

  // On ajoute les colonnes.
  StringUniqueArray columns_name(options()->getNSteps());
  for(Integer i = 0; i < options()->getNSteps(); i++){
    columns_name[i] = "Iteration " + String::fromNumber(i+1);
  }
  m_csv->addColumns(columns_name);

  // On organise les lignes dans le tableau de résultat.
  m_csv->addRows(StringUniqueArray {
    "Sampling duration (Proc)",
    "Tracking duration (Proc)",
    "m_start (Proc)",
    "m_source (Proc)",
    "m_rr (Proc)", 
    "m_split (Proc)",
    "m_absorb (Proc)",
    "m_scatter (Proc)",
    "m_fission (Proc)",
    "m_produce (Proc)",
    "m_collision (Proc)", 
    "m_escape (Proc)",
    "m_census (Proc)",
    "m_num_segments (Proc)", 
    "m_end (Proc)", 
    "m_incoming (Proc)",
    "m_outgoing (Proc)",
    "sum_scalar_flux_tally (Proc)",
  });

  if(options()->getLoadBalancingMat() || options()->getLoadBalancingLoop() > 0){
    m_csv->addRows(StringUniqueArray {
    "LB Sum Criterion Cell Before (Proc)",
    "LB Sum Criterion Cell After (Proc)",
    "LB Sum Criterion Face Before (Proc)",
    "LB Sum Criterion Face After (Proc)",
    });
  }

  if(parallelMng()->commSize() != 1 && parallelMng()->commRank() == 0) {
    m_csv->addRows(StringUniqueArray {
      "Sampling duration (ReduceMax)",
      "Tracking duration (ReduceMax)",
      "m_start (ReduceSum)","m_start (ReduceMin)","m_start (ReduceMax)","m_start (ReduceAvg)",
      "m_source (ReduceSum)","m_source (ReduceMin)","m_source (ReduceMax)","m_source (ReduceAvg)",
      "m_rr (ReduceSum)","m_rr (ReduceMin)","m_rr (ReduceMax)","m_rr (ReduceAvg)",
      "m_split (ReduceSum)","m_split (ReduceMin)","m_split (ReduceMax)","m_split (ReduceAvg)",
      "m_absorb (ReduceSum)","m_absorb (ReduceMin)","m_absorb (ReduceMax)","m_absorb (ReduceAvg)",
      "m_scatter (ReduceSum)","m_scatter (ReduceMin)","m_scatter (ReduceMax)","m_scatter (ReduceAvg)",
      "m_fission (ReduceSum)","m_fission (ReduceMin)","m_fission (ReduceMax)","m_fission (ReduceAvg)",
      "m_produce (ReduceSum)","m_produce (ReduceMin)","m_produce (ReduceMax)","m_produce (ReduceAvg)",
      "m_collision (ReduceSum)","m_collision (ReduceMin)","m_collision (ReduceMax)","m_collision (ReduceAvg)",
      "m_escape (ReduceSum)","m_escape (ReduceMin)","m_escape (ReduceMax)","m_escape (ReduceAvg)",
      "m_census (ReduceSum)","m_census (ReduceMin)","m_census (ReduceMax)","m_census (ReduceAvg)",
      "m_num_segments (ReduceSum)","m_num_segments (ReduceMin)","m_num_segments (ReduceMax)","m_num_segments (ReduceAvg)",
      "m_end (ReduceSum)","m_end (ReduceMin)","m_end (ReduceMax)","m_end (ReduceAvg)",
      "m_incoming (ReduceSum)","m_incoming (ReduceMin)","m_incoming (ReduceMax)","m_incoming (ReduceAvg)",
      "m_outgoing (ReduceSum)","m_outgoing (ReduceMin)","m_outgoing (ReduceMax)","m_outgoing (ReduceAvg)",
      "sum_scalar_flux_tally (ReduceSum)","sum_scalar_flux_tally (ReduceMin)","sum_scalar_flux_tally (ReduceMax)","sum_scalar_flux_tally (ReduceAvg)"
    });
  }
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
 * @brief Méthode permettant d'effectuer l'équilibrage de charge pré-boucle
 * (si l'option loadBalancingMat == true).
 */
void QSModule::
startLoadBalancing()
{
  if(!options()->getLoadBalancingMat()) return;
  info() << "startLoadBalancing";
  loadBalancing();
}

/**
 * @brief Méthode permettant d'effectuer l'équilibrage de charge post-boucle
 * (si itération % option loadBalancingLoop == 0).
 */
void QSModule::
loopLoadBalancing()
{
  if(options()->getLoadBalancingLoop() == 0 
    || (m_global_iteration() % options()->getLoadBalancingLoop() != 0 && m_global_iteration() != 1)
    || m_global_iteration() == options()->getNSteps()) return;
  info() << "loopLoadBalancing";
  loadBalancing();
}

/**
 * @brief Méthode permettant d'effectuer l'équilibrage de charge.
 */
void QSModule::
loadBalancing()
{
  Real sum_cell = 0, sum_face = 0;

  ENUMERATE_CELL(icell, ownCells()){
    sum_cell += m_criterion_lb_cell[icell];
  }

  ENUMERATE_FACE(iface, ownFaces()){
    sum_face += m_criterion_lb_face[iface];
  }

  //parallelMng()->barrier();
  //pinfo() << "P" << mesh()->parallelMng()->commRank() << " - Load Balancing - Difficulté SD avant LB - Cell : " << sum_cell << " - Face : " << sum_face;

  m_csv->editElem(("Iteration " + String::fromNumber(m_global_iteration()+1)), "LB Sum Criterion Cell Before (Proc)", sum_cell);
  m_csv->editElem(("Iteration " + String::fromNumber(m_global_iteration()+1)), "LB Sum Criterion Face Before (Proc)", sum_face);

  ILoadBalanceMng* lb = subDomain()->loadBalanceMng();
  lb->addCriterion(m_criterion_lb_cell);
  lb->addCommCost(m_criterion_lb_face);
  subDomain()->timeLoopMng()->registerActionMeshPartition((IMeshPartitionerBase*)options()->partitioner());
}

/**
 * @brief Méthode permettant d'effectuer l'après équilibrage de charge.
 */
void QSModule::
afterLoadBalancing()
{
  if(
    (options()->getLoadBalancingLoop() == 0 || m_global_iteration() % options()->getLoadBalancingLoop() != 0)
    && !options()->getLoadBalancingMat()
    ) return;

  info() << "AfterLoadBalancing";

  Real sum_cell = 0, sum_face = 0;

  ENUMERATE_CELL(icell, ownCells()){
    sum_cell += m_criterion_lb_cell[icell];
  }

  ENUMERATE_FACE(iface, ownFaces()){
    sum_face += m_criterion_lb_face[iface];
  }

  //pinfo() << "P" << mesh()->parallelMng()->commRank() << " - Load Balancing - Difficulté SD après LB - Cell : " << sum_cell << " - Face : " << sum_face;

  m_csv->editElem(("Iteration " + String::fromNumber(m_global_iteration())), "LB Sum Criterion Cell After (Proc)", sum_cell);
  m_csv->editElem(("Iteration " + String::fromNumber(m_global_iteration())), "LB Sum Criterion Face After (Proc)", sum_face);

  m_criterion_lb_cell.fill(1.);
  m_criterion_lb_face.fill(1.);
}

/**
 * @brief Méthode appelée à la fin de la boucle en temps.
 */
void QSModule::
endModule()
{
  if(parallelMng()->commRank() == 0) {
    RealUniqueArray max_tracking_times(m_csv->getRow("Tracking duration (ReduceMax)"));
    RealUniqueArray num_segments(m_csv->getRow(
      (parallelMng()->commSize() != 1) ?
      "m_num_segments (ReduceSum)" :
      "m_num_segments (Proc)"
    ));

    Real sum_times = 0;
    Int64 sum_segs = 0;

    for(Integer i = 0; i < max_tracking_times.size(); i++) {
      sum_times += max_tracking_times[i];
      sum_segs += num_segments[i];
    }

    Real fOm = sum_segs / sum_times;

    info() << "----------------------------------------";
    info() << " --- Figure Of Merit : "
           << fOm
           << " [Num Segments / Cycle Tracking Time] ---";
    info() << "----------------------------------------";

    m_csv->addElemRow("Figure Of Merit", fOm);
  }
  
  if(options()->getCsvName() != "" || options()->getCsvDir() != "") {
    String path;
    if(options()->getCsvDir() != "")
      path = options()->getCsvDir();
    else
      path = "./csv_output/";

    //m_csv->print();
    m_csv->writeFile(path);
    
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

  m_pre_lb = options()->getLoadBalancingMat();
  m_loop_lb = options()->getLoadBalancingLoop();

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