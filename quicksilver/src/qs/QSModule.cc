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

#include <arcane/utils/ApplicationInfo.h>
#include <arcane/utils/CommandLineArguments.h>

#include <iostream>
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Méthode permettant de lancer l'initialisation des grandeurs
 * au maillage et le service CSV.
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

  // On récupère un pointeur vers le singleton csv.
  m_csv = ServiceBuilder<ISimpleTableOutput>(subDomain()).getSingleton();
  m_csv_compare = ServiceBuilder<ISimpleTableComparator>(subDomain()).getSingleton();

  // Initialisation de la sortie CSV.
  String csvName, csvDir;

  if(options()->getCsvName() != "")
    csvName = options()->getCsvName();
  else
    csvName = "QAMA_P@proc_id@";

  // On regarde si on a le nom du répertoire.
  // Sinon, on donne une emplacement par défaut.
  if(options()->getCsvDir() != "")
    csvDir = options()->getCsvDir();
  else
    csvDir = "csv_output";

  m_csv->init(csvName, csvDir);


  // On ajoute les colonnes (une par itération).
  StringUniqueArray columns_name(options()->getNSteps());
  for(Integer i = 0; i < options()->getNSteps(); i++){
    columns_name[i] = "Iteration " + String::fromNumber(i+1);
  }
  m_csv->addColumns(columns_name);

  // On organise les lignes dans le tableau de résultat.
  // On commence par les infos au niveau du proc.
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

  // Si on a activé l'équilibrage, on ajoute des lignes dédiées.
  if(options()->getLoadBalancingMat() || options()->getLoadBalancingLoop() > 0){
    m_csv->addRows(StringUniqueArray {
      "LB Sum Criterion Cell Before (Proc)",
      "LB Sum Criterion Cell After (Proc)",
      "LB Sum Criterion Face Before (Proc)",
      "LB Sum Criterion Face After (Proc)",
    });
  }

  if(parallelMng()->commRank() == 0) {
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
 * @brief Méthode permettant d'afficher les informations de fin d'itération et de voir
 * si l'on a atteint la dernière itération.
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
 * @brief Méthode permettant d'effectuer l'équilibrage de charge (à la fin de l'itération actuelle).
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

  // Pour les infos d'équilibrage, on affiche l'équilibrage effectuée avant l'itération (il y a donc un décallage ici).
  m_csv->editElement(("Iteration " + String::fromNumber(m_global_iteration()+1)), "LB Sum Criterion Cell Before (Proc)", sum_cell);
  m_csv->editElement(("Iteration " + String::fromNumber(m_global_iteration()+1)), "LB Sum Criterion Face Before (Proc)", sum_face);

  ILoadBalanceMng* lb = subDomain()->loadBalanceMng();
  lb->addCriterion(m_criterion_lb_cell);
  lb->addCommCost(m_criterion_lb_face);
  subDomain()->timeLoopMng()->registerActionMeshPartition((IMeshPartitionerBase*)options()->partitioner());
}

/**
 * @brief Méthode permettant d'effectuer l'après équilibrage de charge (au début de l'itération suivante).
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

  // Pour les infos d'équilibrage, on affiche l'équilibrage effectuée avant le début de l'itération.
  // (ici en revanche, pas de décallages car cette méthode est appelée au début de l'itération d'après).
  m_csv->editElement(("Iteration " + String::fromNumber(m_global_iteration())), "LB Sum Criterion Cell After (Proc)", sum_cell);
  m_csv->editElement(("Iteration " + String::fromNumber(m_global_iteration())), "LB Sum Criterion Face After (Proc)", sum_face);

  m_criterion_lb_cell.fill(1.);
  m_criterion_lb_face.fill(1.);
}

/**
 * @brief Méthode permettant de comparer les résultats obtenus avec
 * les résultats de références.
 */
void QSModule::
compareWithReference()
{
  String reference_input = subDomain()->applicationInfo().commandLineArguments().getParameter("ReferenceDirectory");
  bool overwrite_reference = (subDomain()->applicationInfo().commandLineArguments().getParameter("OverwriteReference") == "true");

  // L'argument overwrite ne fonctionne qu'avec l'argument ref dir.
  if(reference_input.empty()) {
    reference_input = options()->getCsvReferenceDir();
    overwrite_reference = options()->getCsvOverwriteReference();
  }

  if(!reference_input.empty()) {
    info() << "---------------------------------------";
    info() << "-----------Comparator part-------------";
    info() << "---------------------------------------";
    m_csv_compare->init(m_csv);

    if(reference_input != "default") {
      info() << "  Set reference directory: " << reference_input;
      m_csv_compare->editRootDirectory(Directory(reference_input));
    }

    // Si demande d'écriture.
    if(overwrite_reference) {
      info() << "  Write reference file";
      m_csv_compare->writeReferenceFile(0);
    }

    // Sinon lecture.
    else {
      // Si le fichier existe, comparaison.
      if(m_csv_compare->isReferenceExist(0)) {
        m_csv_compare->editRegexRows("^.*ReduceSum.*$");
        m_csv_compare->addRowForComparing("m_incoming (ReduceSum)");
        m_csv_compare->addRowForComparing("m_outgoing (ReduceSum)");
        m_csv_compare->isAnArrayExclusiveRows(true);

        info() << "  Check results with reference file";
        if(!m_csv_compare->compareWithReference(0, 0.01, false)){
          error() << "End checking : Differents values found";
          ARCANE_FATAL("End checking : Differents values found");
        }

        else if(parallelMng()->commRank() == 0){
          info() << "  End checking : Same values!!!";
        }
      }
      // Sinon erreur.
      else {
        error() << "  Reference file not found";
        ARCANE_FATAL("Reference file not found");
      }
    }
    info() << "---------------------------------------";
    info() << "---------End Comparator part-----------";
    info() << "---------------------------------------";
    info() << "-";
  }
}

/**
 * @brief Méthode appelée à la fin de la boucle en temps.
 */
void QSModule::
endModule()
{
  // On calcule la valeur "Figure Of Merit". Le calcul est le même que dans QS originale.
  if(parallelMng()->commRank() == 0) {
    // On utilise les valeurs enregistrées dans le csv.
    RealUniqueArray max_tracking_times(m_csv->row(
      "Tracking duration (ReduceMax)"
    ));
    RealUniqueArray num_segments(m_csv->row(
      "m_num_segments (ReduceSum)"
    ));

    Real sum_times = 0;
    Int64 sum_segs = 0;

    // On fait la somme des temps et des segments (1 segment = 1 tour de tracking).
    for(Integer i = 0; i < max_tracking_times.size(); i++) {
      sum_times += max_tracking_times[i];
      sum_segs += num_segments[i];
    }

    Real fOm = sum_segs / sum_times;

    // On fait un joli affichage.
    info() << "----------------------------------------------";
    info() << "----------------Figure Of Merit---------------";
    info() << "-----[Num Segments / Cycle Tracking Time]-----";
    info() << "----------------------------------------------";
    info() << "--- Num Segments        : " << sum_segs;
    info() << "--- Cycle Tracking Time : " << sum_times;
    info() << "----------------------------------------------";
    info() << "--- Figure Of Merit     : " << fOm;
    info() << "----------------------------------------------";

    // On ajoute une ligne "de séparation" dans le csv, puis le "Figure Of Merit".
    m_csv->addRow("---------------");
    m_csv->addElementInRow("Figure Of Merit", fOm);
  }

  // Si une des options est édité dans le .arc.
  // À noter que le nom du fichier .csv est le nom du tableau (initialisé dans initModule()).
  if(options()->getCsvName() != "" || options()->getCsvDir() != "") {
    info() << "Write results in CSV files";
    if(!m_csv->writeFile()) error() << "Error write CSV";
    info() << "End write results in CSV files";
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

  // Largeur totale du maillage.
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

    m_cell_center_coord[icell] = 0.0;

    ENUMERATE_FACE (iface, cell.faces()) {
      Face face = *iface;

      m_face_center_coord[iface] = 0.0;

      ENUMERATE_NODE (inode, face.nodes()) {
        m_face_center_coord[iface] += node_coord[inode];
      }

      m_face_center_coord[iface] /= 4;
      m_cell_center_coord[icell] += m_face_center_coord[iface];
    }

    m_cell_center_coord[icell] /= 6;

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
    ARCANE_ASSERT(false, ("Boundary condition undefined."));
    return ParticleEvent::undefined;
  }
}
