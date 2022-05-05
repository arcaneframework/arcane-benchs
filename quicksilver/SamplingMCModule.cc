// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "SamplingMCModule.hh"
#include <arcane/Concurrency.h>
#include <map>
#include <set>
#include "MC_RNG_State.hh"
#include "NVTX_Range.hh"
#include "PhysicalConstants.hh"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Méthode permettant d'initialiser les grandeurs aux particules.
 */
void SamplingMCModule::
initModule()
{
  m_particle_family = mesh()->findItemFamily("ArcaneParticles");
  m_particle_family->setHasUniqueIdMap(false);

  m_particle_coord.resize(3);
  m_particle_coord.fill(0.0);

  m_particle_velocity.resize(3);
  m_particle_velocity.fill(0.0);

  m_particle_dir_cos.resize(3);
  m_particle_dir_cos.fill(0.0);

  m_timer = new Timer(subDomain(), "SamplingMC", Timer::TimerReal);
}

/**
 * @brief Méthode permettant de créer/tuer des particules.
 */
void SamplingMCModule::
cycleInit()
{
  {
    Timer::Sentry ts(m_timer);

    clearCrossSectionCache();
    m_processingView = m_particle_family->view();

    arcaneParallelForeach(m_processingView, [&](ParticleVectorView particles) {
      ENUMERATE_PARTICLE (ipartic, particles) {
        if (m_particle_species[ipartic] != ParticleState::exitedParticle 
        && m_particle_species[ipartic] != ParticleState::censusParticle) {
          ARCANE_FATAL("Particule non traitée dans Sampling");
        }
        m_particle_species[ipartic] = ParticleState::oldParticle;
      }
    });

    m_start_a = m_processingView.size();

    // Création des particules.
    sourceParticles();

    // Réduction ou augmentation du nombre de particules.
    populationControl(); // controls particle population

    pinfo(3) << "P" << mesh()->parallelMng()->commRank()
            << " - SourceParticles: " << m_source_a << " particle(s) created.";
    pinfo(3) << "P" << mesh()->parallelMng()->commRank()
            << " - PopulationControl: " << m_rr_a << " particle(s) killed / "
            << m_split_a << " particle(s) created by splitting.";
    Int64 tmpLog = m_rr_a;

    // Roulette sur les particules avec faible poids.
    rouletteLowWeightParticles(); // Delete particles with low statistical weight

    pinfo(3) << "P" << mesh()->parallelMng()->commRank()
            << " - RouletteLowWeightParticles: " << m_rr_a - tmpLog
            << " particle(s) killed.";

    updateTallies();
  }

  info() << "--- P" << mesh()->parallelMng()->commRank() << " - Sampling duration: " << m_timer->lastActivationTime() << "s ---";
}

/**
 * @brief Méthode appelée à la fin de la boucle en temps.
 */
void SamplingMCModule::
endModule() 
{
  delete(m_timer);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Méthode permettant de remettre à zéro m_total.
 */
void SamplingMCModule::
clearCrossSectionCache()
{
  ENUMERATE_CELL (icell, ownCells()) {
    for (Integer i = 0; i < m_n_groups(); i++) {
      m_total[icell][i] = 0.0;
    }
  }
}

/**
 * @brief Méthode permettant de créer des particules.
 * Le nombre de particules créées dépend du volume de la cellule
 * et du taux de production de particule.
 */
void SamplingMCModule::
sourceParticles()
{
  NVTX_Range range("MC_Source_Now");

  Real local_weight_particles = 0;

  // On regarde le nombre de particule que chaque cellule générera.
  ENUMERATE_CELL (icell, ownCells()) {
    Real cell_weight_particles =
    m_volume[icell] * m_source_rate[icell] * m_global_deltat();
    local_weight_particles += cell_weight_particles;
  }

  Real total_weight_particles = 0;

  total_weight_particles = mesh()->parallelMng()->reduce(
  Parallel::ReduceSum, local_weight_particles);

  Int64 num_particles = options()->getNParticles();

  Real source_fraction = 0.1;
  Real source_particle_weight =
  total_weight_particles / (source_fraction * num_particles);

  // Store the source particle weight for later use.
  m_source_particle_weight = source_particle_weight;

  Int64 particle_count = 0;
  // Int64UniqueArray num_particles_cells_decal(ownCells().size()+1);
  // num_particles_cells_decal[0] = 0;

  // On compte le nombre de particule total qui sera créé.
  ENUMERATE_CELL (icell, ownCells()) {
    // if(icell.index() != icell.localId()) ARCANE_FATAL("Trouve une autre
    // solution");
    Real cell_weight_particles =
    m_volume[icell] * m_source_rate[icell] * m_global_deltat();
    Real cell_num_particles_float =
    cell_weight_particles / source_particle_weight;
    particle_count += (Int64)cell_num_particles_float;
    // num_particles_cells_decal[icell.localId()] = particle_count;
  }

  Int64UniqueArray uids(particle_count);
  Int32UniqueArray local_id_cells(particle_count);
  Int32UniqueArray particles_lid(particle_count);
  std::map<Int64, Int64> rng;
  Integer particle_index_g = 0;

  // On gérère les uniqueId et les graines des futures particules.
  // TODO : On a besoin d'un index global si parallélisation.
  ENUMERATE_CELL (icell, ownCells()) {
    Real cell_weight_particles =
    m_volume[icell] * m_source_rate[icell] * m_global_deltat();
    Real cell_num_particles_float =
    cell_weight_particles / source_particle_weight;
    Integer cell_num_particles = (Integer)cell_num_particles_float;

    // for ( Integer particle_index =
    // num_particles_cells_decal[icell.localId()];
    //      particle_index < num_particles_cells_decal[icell.localId()+1];
    //      particle_index++ )
    for (Integer particle_index = 0; particle_index < cell_num_particles;
         particle_index++) {
      Int64 random_number_seed;
      Int64 rns;
      Int64 id;

      random_number_seed = m_source_tally[icell]; // TODO : Atomic
      m_source_tally[icell]++; // TODO : Atomic

      random_number_seed +=
      (*icell).uniqueId().asInt64() * INT64_C(0x0100000000);

      // La graine sera considérée comme un uint64, mais l'uniqueid doit être
      // positif donc on change le signe en mettant le bit de poids fort à 0.
      rns = rngSpawn_Random_Number_Seed(&random_number_seed);

      id = random_number_seed;
      id &= ~(1UL << 63);

      uids[particle_index_g] = id;
      local_id_cells[particle_index_g] = icell.localId();
      rng[id] = rns;

      particle_index_g++;
    }
  }
  m_particle_family->toParticleFamily()->addParticles(uids, local_id_cells,
                                                      particles_lid);
  m_particle_family->endUpdate();

  ParticleVectorView viewSrcP = m_particle_family->view(particles_lid);

  // Les particules sont créées, on les initialise donc.
  arcaneParallelForeach(viewSrcP, [&](ParticleVectorView particles) {
    ENUMERATE_PARTICLE (ipartic, particles) {
      Particle p = (*ipartic);
      initParticle(p, rng[p.uniqueId().asInt64()]);

      generate3DCoordinate(p);
      sampleIsotropic(p);
      m_particle_kin_ene[p] =
      (m_e_max() - m_e_min()) * rngSample(&m_particle_rns[p]) + m_e_min();

      Real speed = getSpeedFromEnergy(p);

      m_particle_velocity[p][MD_DirX] = speed * m_particle_dir_cos[p][MD_DirA];
      m_particle_velocity[p][MD_DirY] = speed * m_particle_dir_cos[p][MD_DirB];
      m_particle_velocity[p][MD_DirZ] = speed * m_particle_dir_cos[p][MD_DirG];

      m_particle_weight[p] = source_particle_weight;

      Real randomNumber = rngSample(&m_particle_rns[p]);
      m_particle_num_mean_free_path[p] = -1.0 * std::log(randomNumber);

      randomNumber = rngSample(&m_particle_rns[p]);
      m_particle_time_census[p] = m_global_deltat() * randomNumber;

      m_source_a++;
    }
  });

  m_processingView = m_particle_family->view();
}

/**
 * @brief Méthode permettant de controler la quantité de particule.
 * Les particules sont splittées ou tuées selon le nombre de particule
 * ciblé.
 */
void SamplingMCModule::
populationControl()
{
  NVTX_Range range("populationControl");

  Int64 targetNumParticles = options()->getNParticles();
  Int64 globalNumParticles = 0;
  Integer localNumParticles = m_processingView.size();

  globalNumParticles =
  mesh()->parallelMng()->reduce(Parallel::ReduceSum, localNumParticles);

  // Soit on augmente la population (>1), soit on l'a diminue (<1), soit on n'y
  // touche pas (=1).
  Real splitRRFactor = (Real)targetNumParticles / (Real)globalNumParticles;

  // March backwards through the vault so killed particles doesn't mess up the
  // indexing
  if (splitRRFactor == 1) {
    return;
  }
  else if (splitRRFactor < 1) {
    Int32UniqueArray supprP;

    arcaneParallelForeach(m_processingView, [&](ParticleVectorView particles) {
      ENUMERATE_PARTICLE (iparticle, particles) {
        Real randomNumber = rngSample(&m_particle_rns[iparticle]);
        if (randomNumber > splitRRFactor) {
          // Kill
          m_rr_a++;
          GlobalMutex::ScopedLock(m_mutex);
          supprP.add(iparticle.localId());
        }
        else {
          // Ici, splitRRFactor < 1 donc on augmente la taille de la
          // particule.
          m_particle_weight[iparticle] /= splitRRFactor;
        }
      }
    });
    m_particle_family->toParticleFamily()->removeParticles(supprP);
    m_particle_family->toParticleFamily()->endUpdate();
  }

  else if (splitRRFactor > 1) {
    Int64UniqueArray addIdP;
    Int64UniqueArray addRns;
    Int32UniqueArray addCellIdP;
    Int32UniqueArray addSrcP;

    arcaneParallelForeach(m_processingView, [&](ParticleVectorView particles) {
      ENUMERATE_PARTICLE (iparticle, particles) {
        Particle particle = (*iparticle);
        Real randomNumber = rngSample(&m_particle_rns[iparticle]);

        // Split
        Integer splitFactor = (Integer)floor(splitRRFactor);
        if (randomNumber > (splitRRFactor - splitFactor)) {
          splitFactor--;
        }

        m_particle_weight[iparticle] /= splitRRFactor;

        GlobalMutex::ScopedLock(m_mutex);
        for (Integer splitFactorIndex = 0; splitFactorIndex < splitFactor;
             splitFactorIndex++) {
          m_split_a++;
          Int64 rns =
          rngSpawn_Random_Number_Seed(&m_particle_rns[iparticle]);
          addRns.add(rns);
          rns &= ~(1UL << 63); // On passe en positif.
          addIdP.add(rns);
          addCellIdP.add(particle.cell().localId());
          addSrcP.add(iparticle.localId());
        }
      }
    });

    Int32UniqueArray particles_lid(addIdP.size());
    m_particle_family->toParticleFamily()->addParticles(addIdP, addCellIdP,
                                                        particles_lid);
    m_particle_family->toParticleFamily()->endUpdate();

    cloneParticles(addSrcP, particles_lid, addRns);
  }

  m_processingView = m_particle_family->view();
}

/**
 * @brief Méthode permettant d'initialiser les grandeurs de la particule p.
 *
 * @param p La particule à initialiser.
 * @param rns La graine à donner à la particule.
 */
void SamplingMCModule::
initParticle(Particle p, Int64 rns)
{
  m_particle_rns[p] = rns;
  m_particle_coord[p][MD_DirX] = 0.0;
  m_particle_coord[p][MD_DirY] = 0.0;
  m_particle_coord[p][MD_DirZ] = 0.0;

  m_particle_velocity[p][MD_DirX] = 0.0;
  m_particle_velocity[p][MD_DirY] = 0.0;
  m_particle_velocity[p][MD_DirZ] = 0.0;

  m_particle_dir_cos[p][MD_DirA] = 0.0;
  m_particle_dir_cos[p][MD_DirB] = 0.0;
  m_particle_dir_cos[p][MD_DirG] = 0.0;

  m_particle_kin_ene[p] = 0.0;
  m_particle_weight[p] = 0.0;
  m_particle_time_census[p] = 0.0;
  m_particle_total_cross_section[p] = 0.0;
  m_particle_age[p] = 0.0;
  m_particle_num_mean_free_path[p] = 0.0;
  m_particle_mean_free_path[p] = 0.0;
  m_particle_seg_path_length[p] = 0.0;
  m_particle_last_event[p] = ParticleEvent::census;
  m_particle_num_coll[p] = 0;
  m_particle_num_seg[p] = 0.0;
  m_particle_species[p] = ParticleState::newParticle;
  m_particle_ene_grp[p] = 0;
  m_particle_face[p] = 0;
  m_particle_facet[p] = 0;
  m_particle_normal_dot[p] = 0.0;
}

/**
 * @brief Méthode permettant de cloner des particules.
 * Copie des propriétés de particules sources vers les particules destinations.
 * La taille des trois tableaux doit être identique.
 *
 * @param idsSrc Tableau des ids des particules sources.
 * @param idsNew Tableau des ids des particules destinations.
 * @param rnsNew Tableau contenant les graines à donner aux particules.
 */
void SamplingMCModule::
cloneParticles(Int32UniqueArray idsSrc,
               Int32UniqueArray idsNew,
               Int64UniqueArray rnsNew)
{
  ParticleVectorView viewSrcP = m_particle_family->view(idsSrc);
  ParticleVectorView viewNewP = m_particle_family->view(idsNew);
  ENUMERATE_PARTICLE (iparticle, viewSrcP) {
    Particle pNew(viewNewP[iparticle.index()].internal());
    cloneParticle((*iparticle), pNew, rnsNew[iparticle.index()]);
  }
}

/**
 * @brief Méthode permettant de cloner une particule.
 * Copie des propriétés d'une particule source vers une particule destination.
 *
 * @param pSrc La particule source.
 * @param pNew La particule destination.
 * @param rns La graine à donner à la particule destination.
 */
void SamplingMCModule::
cloneParticle(Particle pSrc, Particle pNew, Int64 rns)
{
  m_particle_rns[pNew] = rns;

  m_particle_coord[pNew][MD_DirX] = m_particle_coord[pSrc][MD_DirX];
  m_particle_coord[pNew][MD_DirY] = m_particle_coord[pSrc][MD_DirY];
  m_particle_coord[pNew][MD_DirZ] = m_particle_coord[pSrc][MD_DirZ];

  m_particle_velocity[pNew][MD_DirX] = m_particle_velocity[pSrc][MD_DirX];
  m_particle_velocity[pNew][MD_DirY] = m_particle_velocity[pSrc][MD_DirY];
  m_particle_velocity[pNew][MD_DirZ] = m_particle_velocity[pSrc][MD_DirZ];

  m_particle_dir_cos[pNew][MD_DirA] = m_particle_dir_cos[pSrc][MD_DirA];
  m_particle_dir_cos[pNew][MD_DirB] = m_particle_dir_cos[pSrc][MD_DirB];
  m_particle_dir_cos[pNew][MD_DirG] = m_particle_dir_cos[pSrc][MD_DirG];

  m_particle_kin_ene[pNew] = m_particle_kin_ene[pSrc];
  m_particle_weight[pNew] = m_particle_weight[pSrc];
  m_particle_time_census[pNew] = m_particle_time_census[pSrc];
  m_particle_total_cross_section[pNew] = m_particle_total_cross_section[pSrc];
  m_particle_age[pNew] = m_particle_age[pSrc];
  m_particle_num_mean_free_path[pNew] = m_particle_num_mean_free_path[pSrc];
  m_particle_mean_free_path[pNew] = m_particle_mean_free_path[pSrc];
  m_particle_seg_path_length[pNew] = m_particle_seg_path_length[pSrc];
  m_particle_last_event[pNew] = m_particle_last_event[pSrc];
  m_particle_num_coll[pNew] = m_particle_num_coll[pSrc];
  m_particle_num_seg[pNew] = m_particle_num_seg[pSrc];
  m_particle_species[pNew] = ParticleState::clonedParticle;
  m_particle_ene_grp[pNew] = m_particle_ene_grp[pSrc];
  m_particle_face[pNew] = m_particle_face[pSrc];
  m_particle_facet[pNew] = m_particle_facet[pSrc];
  m_particle_normal_dot[pNew] = m_particle_normal_dot[pSrc];
}

/**
 * @brief Méthode permettant de tuer (aléatoirement) les particules trop
 * petites. Roulette russe.
 */
void SamplingMCModule::
rouletteLowWeightParticles()
{
  NVTX_Range range("rouletteLowWeightParticles");

  const Real lowWeightCutoff = options()->getLowWeightCutoff();

  if (lowWeightCutoff > 0.0) {
    Int32UniqueArray supprP;

    // March backwards through the vault so killed particles don't mess up the
    // indexing
    const Real weightCutoff = lowWeightCutoff * m_source_particle_weight;

    arcaneParallelForeach(m_processingView, [&](ParticleVectorView particles) {
      ENUMERATE_PARTICLE (iparticle, particles) {
        if (m_particle_weight[iparticle] <= weightCutoff) {
          Real randomNumber = rngSample(&m_particle_rns[iparticle]);
          if (randomNumber <= lowWeightCutoff) {
            // The particle history continues with an increased weight.
            m_particle_weight[iparticle] /= lowWeightCutoff;
          }
          else {
            // Kill
            m_rr_a++;
            GlobalMutex::ScopedLock(m_mutex);
            supprP.add(iparticle.localId());
          }
        }
      }
    });

    m_particle_family->toParticleFamily()->removeParticles(supprP);
    m_particle_family->toParticleFamily()->endUpdate();

    m_processingView = m_particle_family->view();
  }
}

/**
 * @brief Méthode permettant de générer des coordonnées aléatoire pour une
 * particule.
 *
 * @param p La particule ayant besoin de ces nouvelles coordonnées.
 */
void SamplingMCModule::
generate3DCoordinate(Particle p)
{
  Cell cell = p.cell();
  Int64* random_number_seed = &m_particle_rns[p];

  // Determine the cell-center nodal point coordinates.
  MC_Vector center(m_coord_center[cell]);

  Real random_number = rngSample(random_number_seed);
  Real which_volume = random_number * 6.0 * m_volume[cell];

  // Find the tet to sample from.
  Real current_volume = 0.0;
  Node first_node;
  Node second_node;
  Face face;

  // Pour pouvoir comparer les résultats avec ceux de QS original,
  // on doit explorer les facet de la même manière que QS original.
  for (Integer i = 0; i < 6; i++) {
    face = cell.face(QS_to_arcaneFace[i]);

    for (Integer j = 0; j < 4; j++) {
      Integer first_pos_node = QS_to_arcaneNode[i * 4 + j];
      Integer second_pos_node =
      QS_to_arcaneNode[i * 4 + ((j == 3) ? 0 : j + 1)];

      first_node = face.node(first_pos_node);
      second_node = face.node(second_pos_node);

      MC_Vector point0(m_coord_cm[first_node]);
      MC_Vector point1(m_coord_cm[second_node]);
      MC_Vector point2(m_coord_mid_cm[face]);

      Real subvolume = computeTetVolume(point0, point1, point2, center);
      current_volume += subvolume;

      if (current_volume >= which_volume)
        break;
    }
    if (current_volume >= which_volume)
      break;
  }
  /*
  ENUMERATE_FACE(iface, cell.faces())
  {
    face = (*iface);
    for (Integer i = 0; i < 4; i++)
    {
      first_node  = face.node(i);
      second_node = face.node(((i == 3) ? 0 : i+1));

      MC_Vector point0(m_coord_cm[first_node]);
      MC_Vector point1(m_coord_cm[second_node]);
      MC_Vector point2(m_coord_mid_cm[face]);

      Real subvolume = computeTetVolume(point0, point1, point2, center);
      current_volume += subvolume;

      if(current_volume >= which_volume) break;
    }
    if(current_volume >= which_volume) break;
  }
  */

  // Sample from the tet.
  Real r1 = rngSample(random_number_seed);
  Real r2 = rngSample(random_number_seed);
  Real r3 = rngSample(random_number_seed);

  // Cut and fold cube into prism.
  if (r1 + r2 > 1.0) {
    r1 = 1.0 - r1;
    r2 = 1.0 - r2;
  }
  // Cut and fold prism into tetrahedron.
  if (r2 + r3 > 1.0) {
    Real tmp = r3;
    r3 = 1.0 - r1 - r2;
    r2 = 1.0 - tmp;
  }
  else if (r1 + r2 + r3 > 1.0) {
    Real tmp = r3;
    r3 = r1 + r2 + r3 - 1.0;
    r1 = 1.0 - r2 - tmp;
  }

  // numbers 1-4 are the barycentric coordinates of the random point.
  Real r4 = 1.0 - r1 - r2 - r3;

  MC_Vector point0(m_coord_cm[first_node]);
  MC_Vector point1(m_coord_cm[second_node]);
  MC_Vector point2(m_coord_mid_cm[face]);

  m_particle_coord[p][MD_DirX] =
  (r4 * center.x + r1 * point0.x + r2 * point1.x + r3 * point2.x);
  m_particle_coord[p][MD_DirY] =
  (r4 * center.y + r1 * point0.y + r2 * point1.y + r3 * point2.y);
  m_particle_coord[p][MD_DirZ] =
  (r4 * center.z + r1 * point0.z + r2 * point1.z + r3 * point2.z);
}

/**
 * @brief Méthode permettant de calculer le volume d'un tetrahédre (x6).
 *
 * @param v0_ Le premier point du tétrahédre.
 * @param v1_ Le second point du tétrahédre.
 * @param v2_ Le troisième point du tétrahédre.
 * @param v3  Le quatrième point du tétrahédre.
 * @return Real Le volume du tétrahédre (x6).
 */
Real SamplingMCModule::
computeTetVolume(const MC_Vector& v0_,
                 const MC_Vector& v1_,
                 const MC_Vector& v2_,
                 const MC_Vector& v3)
{
  MC_Vector v0(v0_), v1(v1_), v2(v2_);

  v0.x -= v3.x;
  v0.y -= v3.y;
  v0.z -= v3.z;
  v1.x -= v3.x;
  v1.y -= v3.y;
  v1.z -= v3.z;
  v2.x -= v3.x;
  v2.y -= v3.y;
  v2.z -= v3.z;

  return v0.z * (v1.x * v2.y - v1.y * v2.x) +
  v0.y * (v1.z * v2.x - v1.x * v2.z) +
  v0.x * (v1.y * v2.z - v1.z * v2.y);
}

/**
 * @brief Méthode permettant de donner une orientation aléatoire à la particule.
 *
 * @param p La particule qui aura la nouvelle orientation.
 */
void SamplingMCModule::
sampleIsotropic(Particle p)
{
  m_particle_dir_cos[p][MD_DirG] = 1.0 - 2.0 * rngSample(&m_particle_rns[p]);
  Real sine_gamma = sqrt((
  1.0 - (m_particle_dir_cos[p][MD_DirG] * m_particle_dir_cos[p][MD_DirG])));
  Real phi =
  PhysicalConstants::_pi * (2.0 * rngSample(&m_particle_rns[p]) - 1.0);

  m_particle_dir_cos[p][MD_DirA] = sine_gamma * cos(phi);
  m_particle_dir_cos[p][MD_DirB] = sine_gamma * sin(phi);
}

/**
 * @brief Méthode permettant de récupérer la vitesse de la particule avec son
 * énergie cinétique.
 *
 * @param p La particule concernée.
 * @return Real La vitesse de la particule.
 */
Real SamplingMCModule::
getSpeedFromEnergy(Particle p)
{
  Real energy = m_particle_kin_ene[p];
  static const Real rest_mass_energy =
  PhysicalConstants::_neutronRestMassEnergy;
  static const Real speed_of_light = PhysicalConstants::_speedOfLight;

  return speed_of_light *
  sqrt(energy * (energy + 2.0 * (rest_mass_energy)) /
       ((energy + rest_mass_energy) * (energy + rest_mass_energy)));
}

/**
 * @brief Méthode permettant de copier les atomics dans les variables Arcane.
 *
 */
void SamplingMCModule::
updateTallies()
{
  m_start = m_start_a;
  m_source = m_source_a;
  m_rr = m_rr_a;
  m_split = m_split_a;

  m_start_a = 0;
  m_source_a = 0;
  m_rr_a = 0;
  m_split_a = 0;
}
