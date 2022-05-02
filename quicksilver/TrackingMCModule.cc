#include "TrackingMCModule.hh"
#include <set>
#include <map>
#include "PhysicalConstants.hh"
#include "MC_RNG_State.hh"
#include "MC_Facet_Geometry.hh"
#include "arcane/Concurrency.h"
#include <thread>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Méthode permettant d'initialiser le ParticleExchanger
 * et les materiaux du maillage.
 */
void TrackingMCModule::
initModule()
{
  m_cartesian_mesh = ICartesianMesh::getReference(mesh(), true);
  m_material_mng = IMeshMaterialMng::getReference(mesh());

  m_particle_family = mesh()->findItemFamily("ArcaneParticles");

  IParticleExchanger* pe = options()->particleExchanger();
  pe->initialize(m_particle_family);

  // Configuration des materiaux.
  initNuclearData(); 
}

/**
 * @brief Méthode permettant de lancer le suivi des particules.
 */
void TrackingMCModule::
cycleTracking()
{
  computeCrossSection();
  tracking();
  updateTallies();
}

/**
 * @brief Méthode appelée à la fin de la boucle en temps.
 */
void TrackingMCModule::
endModule()
{

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Méthode permettant d'initialiser les matériaux du maillage.
 */
void TrackingMCModule::
initNuclearData()
{
  info() << "Initialisation des matériaux";

  m_nuclearData = new NuclearData(m_n_groups(), m_e_min(), m_e_max());

  Integer num_cross_sections = options()->cross_section().size();
  Integer num_materials = options()->material().size();
  Integer num_geometries = options()->geometry().size();
  Integer num_isotopes = 0;

  std::map<String, Integer> crossSection;
  UniqueArray<Polynomial> polynomials;

  for( Integer i = 0; i < num_cross_sections; i++ )
  {
    String cs_name = options()->cross_section[i].getName();

    Real aa = options()->cross_section[i].getA();
    Real bb = options()->cross_section[i].getB();
    Real cc = options()->cross_section[i].getC();
    Real dd = options()->cross_section[i].getD();
    Real ee = options()->cross_section[i].getE();

    crossSection.insert(std::make_pair(cs_name, i));
    polynomials.add(Polynomial(aa, bb, cc, dd, ee));
  }


  Int32UniqueArray localIdMat[num_geometries];
  ENUMERATE_CELL(icell, ownCells())
  {
    // On prend le dernier materiau de la liste (si des geométries ont des zones communes) comme dans QS original.
    for( Integer i = num_geometries-1; i >= 0; i-- )
    {
      if(isInGeometry(i, (*icell)))
      {
        localIdMat[i].add(icell.localId());
        break;
      }
    }
  }
  
  
  for( Integer i = 0; i < num_materials; i++ )
  {
    String mat_name = options()->material[i].getName();
    num_isotopes += options()->material[i].getNIsotopes();
    MeshMaterialInfo* newMaterial = m_material_mng->registerMaterialInfo(mat_name);
    MeshEnvironmentBuildInfo ebi1(mat_name);
    ebi1.addMaterial(mat_name);
    m_material_mng->createEnvironment(ebi1);
  }

  m_material_mng->endCreate();
  m_nuclearData->_isotopes.reserve( num_isotopes );

  ConstArrayView<IMeshMaterial*> materials = m_material_mng->materials();
  {
    MeshMaterialModifier modifier(m_material_mng);
    for( Integer i = 0; i < num_geometries; i++ )
    {
      String materialName = options()->geometry[i].getMaterial();

      for (Integer j = 0; j < num_materials; j++)
      {
        // Le nom du materiau est "nomEnvironnement_nomMateriau" et nomEnvironnement = nomMateriau.
        if(materials[j]->name() == materialName+"_"+materialName)
        {
          modifier.addCells(materials[j], localIdMat[i]);
          break;
        }
      }
    }
  }


  for( Integer i = 0; i < num_materials; i++ )
  {
    String mat_name = options()->material[i].getName();
    Real mass = options()->material[i].getMass();
    Integer nIsotopes = options()->material[i].getNIsotopes();
    Integer nReactions = options()->material[i].getNReactions();
    Real sourceRate = options()->material[i].getSourceRate();
    if(mat_name == "defaultMaterial" && sourceRate == 0.0) sourceRate = 1e10;
    String fissionCrossSection = options()->material[i].getFissionCrossSection();
    String scatteringCrossSection = options()->material[i].getScatteringCrossSection();
    String absorptionCrossSection = options()->material[i].getAbsorptionCrossSection();
    Real totalCrossSection = options()->material[i].getTotalCrossSection();
    Real fissionCrossSectionRatio = options()->material[i].getFissionCrossSectionRatio();
    Real scatteringCrossSectionRatio = options()->material[i].getScatteringCrossSectionRatio();
    Real absorptionCrossSectionRatio = options()->material[i].getAbsorptionCrossSectionRatio();

    Real nuBar = options()->cross_section[crossSection.at(fissionCrossSection)].getNuBar();

    ENUMERATE_MATCELL(icell, materials[i])
    {
      m_mass[icell] = mass;
      m_source_rate[icell] = sourceRate;
    }
    m_iso_gid.resize(nIsotopes);
    m_atom_fraction.resize(nIsotopes);


    for (Integer iIso=0; iIso<nIsotopes; ++iIso)
    {
      Integer isotopeGid = m_nuclearData->addIsotope(
        nReactions,
        polynomials[crossSection.at(fissionCrossSection)],
        polynomials[crossSection.at(scatteringCrossSection)],
        polynomials[crossSection.at(absorptionCrossSection)],
        nuBar,
        totalCrossSection,
        fissionCrossSectionRatio,
        scatteringCrossSectionRatio,
        absorptionCrossSectionRatio);
      
      // atom_fraction for each isotope is 1/nIsotopes.  Treats all
      // isotopes as equally prevalent.
      ENUMERATE_MATCELL(icell, materials[i])
      {
        m_iso_gid[icell][iIso] = isotopeGid;
        m_atom_fraction[icell][iIso] = 1.0/nIsotopes;
      }
    }
  }
}

/**
 * @brief Méthode permettant de suivre les particules jusqu'a que leur temps
 * de census soit atteint.
 * Cette méthode effectue une itération et plusieurs sous-itérations.
 */
void TrackingMCModule::
tracking()
{
  bool done = false;

  // Ne sert qu'a debug maintenant.
  //m_local_ids_processed.clear();

  // Toutes les particles de m_particle_family sont à suivre.
  ParticleVectorView m_processingView = m_particle_family->view();
  ParticleVectorView inView;
  Int32UniqueArray extraClone;
  Int32UniqueArray incomingClone;
  Arcane::ParallelLoopOptions optionns;
  optionns.setGrainSize(1000);


  #if LOG
  pinfo() << "P" << mesh()->parallelMng()->commRank() << " - Tracking of " << m_processingView.size() << " particles.";
  #endif

  Integer particle_count = 0; // Initialize count of num_particles processed
  Integer iter = 1;

  
  while (!done)
  {
    arcaneParallelForeach(m_processingView, optionns, [&](ParticleVectorView particles){
      ENUMERATE_PARTICLE(iparticle, particles)
      {
        Particle particle = (*iparticle);
        
        cycleTrackingGuts(particle);
        #if LOG
        if(iparticle.index() % 50000 == 0)
        {
          debug() << "--------";
          pinfo() << "P" << mesh()->parallelMng()->commRank() << " - SubIter #" << iter << " - Number of particles processed : " << iparticle.index() << "/" << m_processingView.size();
          debug() << "  m_local_ids_processed : " << m_census_a << " m_local_ids_extra : " << m_local_ids_extra.size();
          debug() << "--------";
        }
        #endif
      }
    });

    particle_count += m_processingView.size();
    #if LOG
    if(m_processingView.size() != 0) pinfo() << "P" << mesh()->parallelMng()->commRank() << " - SubIter #" << iter << " - Number of particles processed : " << m_processingView.size() << "/" << m_processingView.size();
    #endif
    
    if(mesh()->parallelMng()->commSize() > 1 && inView.size() > 0)
    {
      #if LOG
      pinfo() << "P" << mesh()->parallelMng()->commRank() << " - SubIter #" << iter << " - Computing incoming particles";
      #endif
      arcaneParallelForeach(inView, optionns, [&](ParticleVectorView particles){
        ENUMERATE_PARTICLE(iparticle, particles)
        {
          Particle particle = (*iparticle);
          cycleTrackingGuts(particle);

          #if LOG
          if(iparticle.index() % 50000 == 0)
          {
            pinfo() << "P" << mesh()->parallelMng()->commRank() << " - SubIter #" << iter << " - Number of incoming particles processed : " << iparticle.index() << "/" << inView.size();
          }
          #endif
        }
      });
      particle_count += inView.size();
      #if LOG
      pinfo() << "P" << mesh()->parallelMng()->commRank() << " - SubIter #" << iter << " - Number of incoming particles processed : " << inView.size() << "/" << inView.size();
      #endif
    }

    #if LOG
    debug() << "========";
    pinfo() << "P" << mesh()->parallelMng()->commRank() << " - End SubIter #" << iter  << " - Total number of particles processed : " << particle_count;
    debug() << "  m_local_ids_exit : " << m_local_ids_exit.size() << " m_local_ids_extra_cellId : " << m_local_ids_extra_cellId.size()<< " m_local_ids_extra : " << m_local_ids_extra.size();
    debug() << "  m_local_ids_processed : " << m_census_a << " m_local_ids_exit : " << m_local_ids_exit.size();
    debug() << "========";
    #endif

    
    // On retire les particules qui sont sortie du maillage.
    m_particle_family->toParticleFamily()->removeParticles(m_local_ids_exit);
    // endUpdate fait par collisionEventSuite;
    m_local_ids_exit.clear();

    // On effectue la suite des collisions, si besoin.
    collisionEventSuite();

    if(mesh()->parallelMng()->commSize() > 1)
    {
      //if(m_local_ids_out.size() > 1000 || m_local_ids_extra.size() == 0)
      //{
      // On essaye de recevoir tant que quelqu'un bosse encore.
      do{

        // Echange particles.
        IParticleExchanger* pe = options()->particleExchanger();
        pe->beginNewExchange(m_local_ids_out.size());
        pe->exchangeItems(m_local_ids_out.size(), m_local_ids_out, m_rank_out, &m_local_ids_in, 0);

        done = ((mesh()->parallelMng()->reduce(
          Parallel::ReduceMax, 
          m_local_ids_in.size() + m_local_ids_extra.size()
          )) == 0
        );

        #if LOG
        debug() << "/////////";
        debug() << "Proc#" << mesh()->parallelMng()->commRank();
        debug() << "Nb Particles out : " << m_local_ids_out.size();
        debug() << "Nb Particles in : " << m_local_ids_in.size();
        debug() << "/////////";
        #endif

        m_rank_out.clear();
        m_local_ids_out.clear();
      } while(m_local_ids_in.size() == 0 && m_local_ids_extra.size() == 0 && !done);

      incomingClone = m_local_ids_in.clone();
      m_local_ids_in.clear();
      inView = m_particle_family->view(incomingClone);
      //}
    }

    else if(m_local_ids_extra.size() == 0)
    {
      done = true;
    }

    extraClone = m_local_ids_extra.clone();
    m_local_ids_extra.clear();
    
    m_processingView = m_particle_family->view(extraClone);
    iter++;
  }
  m_end_a = m_particle_family->view().size();
}

/**
 * @brief Méthode permettant de copier les atomics dans les variables Arcane.
 */
void TrackingMCModule::
updateTallies()
{
  m_absorb = m_absorb_a;
  m_census = m_census_a;
  m_escape = m_escape_a;
  m_collision = m_collision_a;
  m_fission = m_fission_a;
  m_produce = m_produce_a;
  m_scatter = m_scatter_a;
  m_num_segments = m_num_segments_a;
  m_end = m_end_a;
  
  m_absorb_a = 0;
  m_census_a = 0;
  m_escape_a = 0;
  m_collision_a = 0;
  m_fission_a = 0;
  m_produce_a = 0;
  m_scatter_a = 0;
  m_num_segments_a = 0;
  m_end_a = 0;
}

// Returns true if the specified coordinate in inside the specified
// geometry.  False otherwise

/**
 * @brief Méthode permettant de savoir si une cellule est dans un matériau.
 * On a les dimensions du matériau (en cm) et la position du centre de la cellule
 * (en cm), cette méthode permet de savoir si le centre de la cellule est dans le 
 * matériau.
 * 
 * @param pos La position de la géométrie du matériau dans le fichier d'entrée.
 * @param cell La cellule étudiée.
 * @return true Si la cellule à son centre dans le matériau.
 * @return false Si la cellule n'a pas son centre dans le matériau.
 */
bool TrackingMCModule::
isInGeometry(Integer pos, Cell cell)
{
  bool inside = false;
  switch (options()->geometry[pos].getShape())
  {
    case eShape::BRICK:
      {
        Real xMin = ((options()->geometry[pos].getXMin() == -1.0) ? 0.0 : options()->geometry[pos].getXMin());
        Real xMax = ((options()->geometry[pos].getXMax() == -1.0) ? m_lx() : options()->geometry[pos].getXMax());

        Real yMin = ((options()->geometry[pos].getYMin() == -1.0) ? 0.0 : options()->geometry[pos].getYMin());
        Real yMax = ((options()->geometry[pos].getYMax() == -1.0) ? m_ly() : options()->geometry[pos].getYMax());

        Real zMin = ((options()->geometry[pos].getZMin() == -1.0) ? 0.0 : options()->geometry[pos].getZMin());
        Real zMax = ((options()->geometry[pos].getZMax() == -1.0) ? m_lz() : options()->geometry[pos].getZMax());

        if ((m_coord_center[cell][MD_DirX] >= xMin && m_coord_center[cell][MD_DirX] <= xMax) &&
            (m_coord_center[cell][MD_DirY] >= yMin && m_coord_center[cell][MD_DirY] <= yMax) &&
            (m_coord_center[cell][MD_DirZ] >= zMin && m_coord_center[cell][MD_DirZ] <= zMax) )
          inside = true;
      }
      break;

    case eShape::SPHERE:
      {
        Real xCenter = ((options()->geometry[pos].getXCenter() == -1.0) ? (m_lx()/2) : options()->geometry[pos].getXCenter());
        Real yCenter = ((options()->geometry[pos].getYCenter() == -1.0) ? (m_ly()/2) : options()->geometry[pos].getYCenter());
        Real zCenter = ((options()->geometry[pos].getZCenter() == -1.0) ? (m_lz()/2) : options()->geometry[pos].getZCenter());
        Real radius  = ((options()->geometry[pos].getRadius()  == -1.0) ? (m_lx()/2) : options()->geometry[pos].getRadius() );

        MC_Vector center(xCenter, yCenter, zCenter);
        MC_Vector rr(m_coord_center[cell]);
        if ( (rr-center).Length() <= radius)
          inside = true;
      }

      break;
    default:
      qs_assert(false);
  }
  return inside;
}

/**
 * @brief Méthode permettant d'initialiser des grandeurs de la particule puis de lancer son tracking.
 * 
 * @param particle La particule à suivre.
 */
void TrackingMCModule::
cycleTrackingGuts( Particle particle )
{
  if ( m_particle_time_census[particle] <= 0.0 )
  {
      m_particle_time_census[particle] += m_global_deltat();
  }

  if (m_particle_age[particle] < 0.0) 
  { 
    m_particle_age[particle] = 0.0;
  }

  m_particle_ene_grp[particle] = m_nuclearData->getEnergyGroup(m_particle_kin_ene[particle]);

  // loop over this particle until we cannot do anything more with it on this processor
  cycleTrackingFunction(particle);

  //Make sure this particle is marked as completed
  m_particle_species[particle] = -1;
}

/**
 * @brief Méthode permettant de suivre la particule p.
 * Suivi de la particule jusqu'a que son temps de suivi (census) soit fini ou
 * jusqu'a qu'elle sorte du maillage ou
 * jusqu'a qu'elle se split (dans ce cas, elle sera resuivi à la prochaine sous-itération).
 * 
 * @param particle La particule à suivre.
 */
void TrackingMCModule::
cycleTrackingFunction(Particle particle)
{
  bool keepTrackingThisParticle = false;
  do
  {
      // Determine the outcome of a particle at the end of this segment such as:
      //
      //   (0) Undergo a collision within the current cell,
      //   (1) Cross a facet of the current cell,
      //   (2) Reach the end of the time step and enter census,
      //
      computeNextEvent(particle);
      m_num_segments_a++;

      m_particle_num_seg[particle] += 1.;  /* Track the number of segments this particle has
                                          undergone this cycle on all processes. */
      switch (m_particle_last_event[particle]) {
      case ParticleEvent::collision:
        {

          switch (collisionEvent(particle))
          {
          case 0: // La particule est absorbée.
            {
              GlobalMutex::ScopedLock(m_mutex_exit);
              m_local_ids_exit.add(particle.localId());
            }
            keepTrackingThisParticle = false;
            break;

          case 1: // La particule a juste changée de trajectoire.
            keepTrackingThisParticle = true;
            break;
          
          default: // La particule splitte.
            // On arrete pour pouvoir cloner la particle source.
            keepTrackingThisParticle = false;
            break;
          }
        }
        break;
  
      // On a la particule sur la face, plus qu'à déterminer la suite.
      case ParticleEvent::faceEventUndefined:
        {
          facetCrossingEvent(particle);

          switch (m_particle_last_event[particle])
          {
          case ParticleEvent::cellChange:
            keepTrackingThisParticle = true;
            break;

          case ParticleEvent::escape:
            m_escape_a++;
            m_particle_species[particle] = -1;
            keepTrackingThisParticle = false;
            {
              GlobalMutex::ScopedLock(m_mutex_exit);
              m_local_ids_exit.add(particle.localId());
            }
            break;

          case ParticleEvent::reflection:
            reflectParticle(particle);
            keepTrackingThisParticle = true;
            break;
          
          default:
            // Enters an adjacent cell in an off-processor domain.
            keepTrackingThisParticle = false;
            // Pas de m_local_ids_exit car la particle sera retirée de la famille par ExchangeParticles.
            break;
          }
        }
        break;
  
      case ParticleEvent::census:
        {
          // The particle has reached the end of the time step.
          // {
          //   GlobalMutex::ScopedLock(m_mutex_processed);
          //   m_local_ids_processed.add(particle.localId());
          // }
          m_census_a++;
          keepTrackingThisParticle = false;
          break;
        }
          
      default:
        qs_assert(false);
        break;  // should this be an error
      }
  
  } while ( keepTrackingThisParticle );
}

/**
 * @brief Méthode permettant de finir une collision.
 * Une collision peut aboutir en un dédoublement de particule.
 * Dans ce cas, la particule doit être dédoublé en dehors de l'ENUMERATE_PARTICLE
 * pour que la vue ne soit pas modifié 'en vol' (et donc que le ENUMERATE reste valide).
 */
void TrackingMCModule::
collisionEventSuite()
{
  // On créé les particules.
  Int32UniqueArray particles_lid(m_local_ids_extra_gId.size());

  m_particle_family->toParticleFamily()->addParticles(m_local_ids_extra_gId, m_local_ids_extra_cellId, particles_lid);
  m_particle_family->toParticleFamily()->endUpdate();

  // On leur donne les bonnes propriétés.
  cloneParticles(m_local_ids_extra_srcP, particles_lid, m_local_ids_extra_rns);

  ParticleVectorView viewP = m_particle_family->view(particles_lid);

  // On effectue la suite de la collision.
  for(Integer i = 0; i < particles_lid.size(); i++)
  {
    Particle particle = Particle(viewP[i].internal());
    updateTrajectory(m_local_ids_extra_energyOut[i], m_local_ids_extra_angleOut[i], particle);
  }

  viewP = m_particle_family->view(m_local_ids_extra);

  for(Integer i = 0; i < m_local_ids_extra.size(); i++)
  {
    Particle particle = Particle(viewP[i].internal());
    updateTrajectory(m_local_ids_extra_energyOut_pSrc[i], m_local_ids_extra_angleOut_pSrc[i], particle);
    m_particle_ene_grp[particle] = m_nuclearData->getEnergyGroup(m_particle_kin_ene[particle]);
  }

  //m_local_ids_extra.copy(particles_lid);

  //m_local_ids_extra.resize(m_local_ids_extra.size() + particles_lid.size());
  // On fusionne la liste des particules sources et la liste des nouvelles particules.
  for(Integer i = 0; i < particles_lid.size(); i++)
  {
    m_local_ids_extra.add(particles_lid[i]);
  }

  m_local_ids_extra_rns.clear();
  m_local_ids_extra_gId.clear();
  m_local_ids_extra_cellId.clear();
  m_local_ids_extra_srcP.clear();
  m_local_ids_extra_energyOut.clear();
  m_local_ids_extra_angleOut.clear();
  m_local_ids_extra_energyOut_pSrc.clear();
  m_local_ids_extra_angleOut_pSrc.clear();
}

/**
 * @brief Méthode permettant de trouver le prochain événement de la particule.
 * On a le choix entre Census, Collision et faceEvent.
 * Census : La particule a fini son itération.
 * Collision : La particule fait une collision avec son environnement.
 * faceEvent : La particule se trouve sur le bord de la cellule et doit faire une action.
 * On choisi l'événement le plus proche.
 * 
 * @param particle La particule à étudier.
 */
void TrackingMCModule::
computeNextEvent(Particle particle)
{
  // initialize distances to large number
  Integer number_of_events = 3;
  RealUniqueArray distance(3);
  distance[0] = distance[1] = distance[2] = 1e80;

  // Calculate the particle speed
  MC_Vector velo(m_particle_velocity[particle]);
  Real particle_speed = velo.Length();

  // Force collision if a census event narrowly preempts a collision
  Integer force_collision = 0 ;
  if ( m_particle_num_mean_free_path[particle] < 0.0 )
  {
    force_collision = 1 ;

    if ( m_particle_num_mean_free_path[particle] > -900.0 )
    {
      printf(" computeNextEvent: m_particle_num_mean_free_path[particle] > -900.0 \n");
    }

    m_particle_num_mean_free_path[particle] = PhysicalConstants::_smallDouble;
  }

  // Randomly determine the distance to the next collision
  // based upon the composition of the current cell.
  //Real macroscopic_total_cross_section = weightedMacroscopicCrossSection(particle.cell(), m_particle_ene_grp[particle]);
  Real macroscopic_total_cross_section = m_total[particle.cell()][m_particle_ene_grp[particle]];
    
  // Cache the cross section
  m_particle_total_cross_section[particle] = macroscopic_total_cross_section;
  if (macroscopic_total_cross_section == 0.0)
  {
    m_particle_mean_free_path[particle] = PhysicalConstants::_hugeDouble;
  }
  else
  {
    m_particle_mean_free_path[particle] = 1.0 / macroscopic_total_cross_section;
  }

  if ( m_particle_num_mean_free_path[particle] == 0.0)
  {
    // Sample the number of mean-free-paths remaining before
    // the next collision from an exponential distribution.
    Real random_number = rngSample(&m_particle_rns[particle]);

    m_particle_num_mean_free_path[particle] = -1.0*std::log(random_number);
  }

  // Calculate the distances to collision, nearest facet, and census.

  // Forced collisions do not need to move far.
  if (force_collision)
  {
    distance[ParticleEvent::collision] = PhysicalConstants::_smallDouble;
  }
  else
  {
    distance[ParticleEvent::collision] = m_particle_num_mean_free_path[particle] * m_particle_mean_free_path[particle];
  }

  // process census
  distance[ParticleEvent::census] = particle_speed*m_particle_time_census[particle];


  //  DEBUG  Turn off threshold for now
  //Real distance_threshold = 10.0 * PhysicalConstants::_hugeDouble;
  // Get the current winning distance.
  //Real current_best_distance = PhysicalConstants::_hugeDouble;


  //bool new_segment =  (m_particle_num_seg[particle] == 0 ||
  //                     m_particle_last_event[particle] == ParticleEvent::collision);

  // Calculate the minimum distance to each facet of the cell.
  NearestFacet nearest_facet = getNearestFacet(particle);

  m_particle_normal_dot[particle] = nearest_facet.dot_product;

  distance[ParticleEvent::faceEventUndefined] = nearest_facet.distance_to_facet;

  // Get out of here if the tracker failed to bound this particle's volume.
  if (m_particle_last_event[particle] == ParticleEvent::faceEventUndefined)
  {
    return;
  }

  // Calculate the minimum distance to the selected events.

  // Force a collision (if required).
  if ( force_collision == 1 )
  {
    distance[ParticleEvent::faceEventUndefined] = PhysicalConstants::_hugeDouble;
    distance[ParticleEvent::census]             = PhysicalConstants::_hugeDouble;
    distance[ParticleEvent::collision]          = PhysicalConstants::_tinyDouble;
  }

  // we choose our segment outcome here
  ParticleEvent segment_outcome = (ParticleEvent) findMin(distance);
  

  if (distance[segment_outcome] < 0)
  {
    qs_assert(false);
  }
  
  m_particle_seg_path_length[particle] = distance[segment_outcome];
  m_particle_num_mean_free_path[particle] -= m_particle_seg_path_length[particle] / m_particle_mean_free_path[particle];
  m_particle_last_event[particle] = segment_outcome;

  // Set the segment path length to be the minimum of
  //   (i)   the distance to collision in the cell, or
  //   (ii)  the minimum distance to a facet of the cell, or
  //   (iii) the distance to census at the end of the time step
  if (segment_outcome == ParticleEvent::collision)
  {
    m_particle_num_mean_free_path[particle] = 0.0;
  }

  else if (segment_outcome == ParticleEvent::faceEventUndefined)
  {
    m_particle_face[particle] = nearest_facet.facet / 4;
    m_particle_facet[particle] = nearest_facet.facet % 4;
  }
  
  else if (segment_outcome == ParticleEvent::census)
  {
    m_particle_time_census[particle] = std::min(m_particle_time_census[particle], 0.0);
  }

  // If collision was forced, set mc_particle.num_mean_free_paths = 0
  // so that a new value is randomly selected on next pass.
  if (force_collision == 1)
  {
    m_particle_num_mean_free_path[particle] = 0.0;
  }

  // Do not perform any tallies if the segment path length is zero.
  //   This only introduces roundoff errors.
  if (m_particle_seg_path_length[particle] == 0.0)
  {
    return;
  }

  // Move particle to end of segment, accounting for some physics processes along the segment.

  // Project the particle trajectory along the segment path length.

  m_particle_coord[particle][MD_DirX] += (m_particle_dir_cos[particle][MD_DirA] * m_particle_seg_path_length[particle]);
  m_particle_coord[particle][MD_DirY] += (m_particle_dir_cos[particle][MD_DirB] * m_particle_seg_path_length[particle]);
  m_particle_coord[particle][MD_DirZ] += (m_particle_dir_cos[particle][MD_DirG] * m_particle_seg_path_length[particle]);

  Real segment_path_time = (m_particle_seg_path_length[particle]/particle_speed);

  // Decrement the time to census and increment age.
  m_particle_time_census[particle] -= segment_path_time;
  m_particle_age[particle] += segment_path_time;

  // Ensure mc_particle.time_to_census is non-negative.
  if (m_particle_time_census[particle] < 0.0)
  {
    m_particle_time_census[particle] = 0.0;
  }

  // Accumulate the particle's contribution to the scalar flux.
  // Atomic
  GlobalMutex::ScopedLock(m_mutex_total);
  m_scalar_flux_tally[particle.cell()][m_particle_ene_grp[particle]] += m_particle_seg_path_length[particle] * m_particle_weight[particle];
}

/**
 * @brief Méthode permettant d'exécuter l'événement Collision sur la particule p.
 * 
 * @param particle La particule à modifier.
 * @return Integer 0 si la particule est absorbé par la maille.
 * 1 si la particule fait simplement une collision.
 * 2 ou plus si la particule doit se splitter.
 */
Integer TrackingMCModule::
collisionEvent(Particle particle)
{
  Cell cell = particle.cell();

  //------------------------------------------------------------------------------------------------------------------
  //    Pick the isotope and reaction.
  //------------------------------------------------------------------------------------------------------------------
  Real randomNumber = rngSample(&m_particle_rns[particle]);
  Real totalCrossSection = m_particle_total_cross_section[particle];
  Real currentCrossSection = totalCrossSection * randomNumber;
  Integer selectedIso = -1;
  Integer selectedUniqueNumber = -1;
  Integer selectedReact = -1;
  Integer numIsos = m_iso_gid[cell].size();
  
  for (Integer isoIndex = 0; isoIndex < numIsos && currentCrossSection >= 0; isoIndex++)
  {
    Integer uniqueNumber = m_iso_gid[cell][isoIndex];
    Integer numReacts = m_nuclearData->getNumberReactions(uniqueNumber);
    for (Integer reactIndex = 0; reactIndex < numReacts; reactIndex++)
    {
      currentCrossSection -= macroscopicCrossSection(reactIndex, cell,
                isoIndex, m_particle_ene_grp[particle]);
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
  Real energyOut[MAX_PRODUCTION_SIZE];
  Real angleOut[MAX_PRODUCTION_SIZE];
  Integer nOut = 0;
  Real mat_mass = m_mass[cell];

  m_nuclearData->_isotopes[selectedUniqueNumber]._species[0]._reactions[selectedReact].sampleCollision(
    m_particle_kin_ene[particle], mat_mass, &energyOut[0], &angleOut[0], nOut, &(m_particle_rns[particle]), MAX_PRODUCTION_SIZE );

  //--------------------------------------------------------------------------------------------------------------
  //  Post-Collision Phase 1:
  //    Tally the collision
  //--------------------------------------------------------------------------------------------------------------

  // Set the reaction for this particle.
  m_collision_a++;
  NuclearDataReaction::Enum reactionType = m_nuclearData->_isotopes[selectedUniqueNumber]._species[0].\
          _reactions[selectedReact]._reactionType;

  switch (reactionType)
  {
    case NuclearDataReaction::Scatter:
        m_scatter_a++;
        break;
    case NuclearDataReaction::Absorption:
        m_absorb_a++;
        break;
    case NuclearDataReaction::Fission:
        m_fission_a++;
        m_produce_a += nOut;
        break;
    case NuclearDataReaction::Undefined:
        printf("reactionType invalid\n");
        qs_assert(false);
  }

  if( nOut == 0 ) 
  {
    return 0;
  }

  else if (nOut == 1)
  {
    updateTrajectory(energyOut[0], angleOut[0], particle);
    m_particle_ene_grp[particle] = m_nuclearData->getEnergyGroup(m_particle_kin_ene[particle]);
  }

  else
  {
    GlobalMutex::ScopedLock(m_mutex_extra);
    for (Integer secondaryIndex = 1; secondaryIndex < nOut; secondaryIndex++)
    {
      Int64 rns = rngSpawn_Random_Number_Seed(&m_particle_rns[particle]);
      m_local_ids_extra_rns.add(rns);
      rns &= ~(1UL << 63);
      m_local_ids_extra_gId.add(rns);
      m_local_ids_extra_cellId.add(particle.cell().localId());
      m_local_ids_extra_srcP.add(particle.localId());
      m_local_ids_extra_energyOut.add(energyOut[secondaryIndex]);
      m_local_ids_extra_angleOut.add(angleOut[secondaryIndex]);
    }

    m_local_ids_extra.add(particle.localId());
    m_local_ids_extra_energyOut_pSrc.add(energyOut[0]);
    m_local_ids_extra_angleOut_pSrc.add(angleOut[0]);
  }

  return nOut;
}


/**
 * @brief Méthode permettant d'exécuter l'événement faceEventUndefined.
 * Et de définir le 'Undefined' :
 * Soit escape (sortie de maillage).
 * Soit reflection (rebondie sur la face)
 * Soit cellChange (changement de maille interne au sous-domaine)
 * Soit subDChange (changement de maille externe au sous-domaine)
 * 
 * @param particle La particule à étudier.
 */
void TrackingMCModule::
facetCrossingEvent(Particle particle)
{
  GlobalMutex::ScopedLock(m_mutex_out);
  Face face = particle.cell().face(m_particle_face[particle]);
  m_particle_last_event[particle] = m_boundary_cond[face];

  if ( m_boundary_cond[face] == ParticleEvent::cellChange )
  {
    // The particle will enter into an adjacent cell.
    Cell cell = face.frontCell();

    if(cell == particle.cell())
    {
      cell = face.backCell();
    }
    
    m_particle_family->toParticleFamily()->setParticleCell(particle, cell);
  }
  else if ( m_boundary_cond[face] == ParticleEvent::subDChange )
  {
    // The particle will enter into an adjacent cell on a spatial neighbor.
    
    Cell cell = face.frontCell();

    if(cell == particle.cell())
    {
      cell = face.backCell();
    }
    
    m_particle_family->toParticleFamily()->setParticleCell(particle, cell);

    m_local_ids_out.add(particle.localId());
    m_rank_out.add(cell.owner());
  }
}

/**
 * @brief Méthode permettant d'exécuter l'événement reflection.
 * 
 * @param particle La particule à étudier.
 */
void TrackingMCModule::
reflectParticle(Particle particle)
{
  Integer facet = m_particle_facet[particle];

  Cell cell = particle.cell();
  Face face = cell.face(m_particle_face[particle]);

  Integer first_pos_node = (scan_order[m_index_arc[face]] ? ((facet == 3) ? 0 : facet+1) : facet);
  Integer second_pos_node = (scan_order[m_index_arc[face]] ? facet : ((facet == 3) ? 0 : facet+1));

  Node first_node = face.node(first_pos_node);
  Node second_node = face.node(second_pos_node);

  MC_Vector point0(m_coord_cm[first_node]);
  MC_Vector point1(m_coord_cm[second_node]);
  MC_Vector point2(m_coord_mid_cm[face]);

  MC_General_Plane plane(point0, point1, point2);

  MC_Vector facet_normal(plane.A, plane.B, plane.C);


  Real dot = 2.0*( m_particle_dir_cos[particle][MD_DirA] * facet_normal.x +
                      m_particle_dir_cos[particle][MD_DirB] * facet_normal.y +
                      m_particle_dir_cos[particle][MD_DirG] * facet_normal.z );

  if ( dot > 0 ) // do not reflect a particle that is ALREADY pointing inward
  {
      // reflect the particle
      m_particle_dir_cos[particle][MD_DirA] -= dot * facet_normal.x;
      m_particle_dir_cos[particle][MD_DirB] -= dot * facet_normal.y;
      m_particle_dir_cos[particle][MD_DirG] -= dot * facet_normal.z;
  }

  // Calculate the reflected, velocity components.
  MC_Vector velo(m_particle_velocity[particle]);
  Real particle_speed = velo.Length();
  m_particle_velocity[particle][MD_DirX] = particle_speed * m_particle_dir_cos[particle][MD_DirA];
  m_particle_velocity[particle][MD_DirY] = particle_speed * m_particle_dir_cos[particle][MD_DirB];
  m_particle_velocity[particle][MD_DirZ] = particle_speed * m_particle_dir_cos[particle][MD_DirG];
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
void TrackingMCModule::
cloneParticles(Int32UniqueArray idsSrc, Int32UniqueArray idsNew, Int64UniqueArray rnsNew)
{
  ParticleVectorView viewSrcP = m_particle_family->view(idsSrc);
  ParticleVectorView viewNewP = m_particle_family->view(idsNew);
  ENUMERATE_PARTICLE(iparticle, viewSrcP)
  {
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
void TrackingMCModule::
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
  m_particle_species[pNew] = m_particle_species[pSrc];
  m_particle_ene_grp[pNew] = m_particle_ene_grp[pSrc];
  m_particle_face[pNew] = m_particle_face[pSrc];
  m_particle_facet[pNew] = m_particle_facet[pSrc];
  m_particle_normal_dot[pNew] = m_particle_normal_dot[pSrc];
}

/**
 * @brief Méthode permettant de mettre à jour la trajectoire de la particule.
 * 
 * @param energy Energie cinétique de la particule.
 * @param angle Angle permettant de mettre à jour la trajectoire de la particule.
 * @param particle La particule à traiter.
 */
void TrackingMCModule::
updateTrajectory( Real energy, Real angle, Particle particle )
{
  m_particle_kin_ene[particle] = energy;
  Real cosTheta = angle;
  Real randomNumber = rngSample(&m_particle_rns[particle]);
  Real phi = 2 * 3.14159265 * randomNumber;
  Real sinPhi = sin(phi);
  Real cosPhi = cos(phi);
  Real sinTheta = sqrt((1.0 - (cosTheta*cosTheta)));

  rotate3DVector(particle, sinTheta, cosTheta, sinPhi, cosPhi);

  Real speed = (PhysicalConstants::_speedOfLight *
          sqrt((1.0 - ((PhysicalConstants::_neutronRestMassEnergy *
          PhysicalConstants::_neutronRestMassEnergy) /
          ((energy + PhysicalConstants::_neutronRestMassEnergy) *
          (energy + PhysicalConstants::_neutronRestMassEnergy))))));

  m_particle_velocity[particle][MD_DirX] = speed * m_particle_dir_cos[particle][MD_DirA];
  m_particle_velocity[particle][MD_DirY] = speed * m_particle_dir_cos[particle][MD_DirB];
  m_particle_velocity[particle][MD_DirZ] = speed * m_particle_dir_cos[particle][MD_DirG];

  randomNumber = rngSample(&m_particle_rns[particle]);
  m_particle_num_mean_free_path[particle] = -1.0*std::log(randomNumber);
}

/**
 * @brief Méthode permettant de déterminer les cross sections de toutes les mailles.
 */
void TrackingMCModule::
computeCrossSection()
{
  Arcane::ParallelLoopOptions options;
  //options.setGrainSize(50);

  arcaneParallelForeach(ownCells(), options, [&](CellVectorView cells){
    ENUMERATE_CELL(icell, cells)
    {
      for(Integer i = 0; i < m_n_groups(); i++)
      {
        weightedMacroscopicCrossSection((*icell), i);
      }
    }
  });
}

/**
 * @brief Méthode permettant de déterminer la distance entre la particule p et la prochaine collision.
 * 
 * @param cell La cellule où se trouve la particule
 * @param energyGroup Le groupe d'energie.
 * @return Real La distance entre la particule et la prochaine collision.
 */
void TrackingMCModule::
weightedMacroscopicCrossSection(Cell cell, Integer energyGroup)
{
  // GlobalMutex::ScopedLock(m_mutex_total);
  // Real precomputedCrossSection = m_total[cell][energyGroup];

  // if (precomputedCrossSection > 0.0)
  //   return precomputedCrossSection;
  
  Integer nIsotopes = m_iso_gid[cell].size();
  Real sum = 0.0;
  for (Integer isoIndex = 0; isoIndex < nIsotopes; isoIndex++)
  {
    sum += macroscopicCrossSection(-1, cell, isoIndex, energyGroup);
  }

  m_total[cell][energyGroup] = sum; // Atomic

  //return sum;
}

/**
 * @brief Routine MacroscopicCrossSection calculates the number-density-weighted 
 * macroscopic cross section of a cell.
 * 
 * @param reactionIndex 
 * @param cell 
 * @param isoIndex 
 * @param energyGroup 
 * @return Real 
 */
Real TrackingMCModule::
macroscopicCrossSection(Integer reactionIndex, Cell cell, Integer isoIndex, Integer energyGroup)
{
  // Initialize various data items.

  Real atom_fraction = m_atom_fraction[cell][isoIndex];

  Real microscopicCrossSection = 0.0;
  // The cell number density is the fraction of the atoms in cell
  // volume of this isotope.  We set this (elsewhere) to 1/nIsotopes.
  // This is a statement that we treat materials as if all of their
  // isotopes are present in equal amounts
  Real cell_number_density = m_cell_number_density[cell];

  Integer isotopeGid = m_iso_gid[cell][isoIndex];
  if ( atom_fraction == 0.0 || cell_number_density == 0.0) 
  { 
    return 1e-20;
  }

  if (reactionIndex < 0)
  {
    // Return total cross section
    microscopicCrossSection = m_nuclearData->getTotalCrossSection(isotopeGid, energyGroup);
  }
  else
  {
    // Return the reaction cross section
    microscopicCrossSection = m_nuclearData->getReactionCrossSection(reactionIndex,
              isotopeGid, energyGroup);
  }

  return atom_fraction * cell_number_density * microscopicCrossSection;
}

/**
 * @brief Méthode permettant de trouver la facet la plus proche de la particule p.
 * 
 * @param particle La particule à étudier.
 * @return NearestFacet Les caractéristiques de la facet la plus proche.
 */
NearestFacet TrackingMCModule::
getNearestFacet(Particle particle)
{
  Cell cell = particle.cell();
  MC_Vector *facet_coords[3];
  Integer iteration = 0;
  Real move_factor = 0.5 * PhysicalConstants::_smallDouble;
  NearestFacet nearest_facet;
  Integer retry = 1;

  while (retry) // will break out when distance is found
  {
    // Determine the distance to each facet of the cell.
    // (1e-8 * Radius)^2
    Real plane_tolerance = 1e-16*(m_particle_coord[particle][MD_DirX]*m_particle_coord[particle][MD_DirX] +
                                    m_particle_coord[particle][MD_DirY]*m_particle_coord[particle][MD_DirY] +
                                    m_particle_coord[particle][MD_DirZ]*m_particle_coord[particle][MD_DirZ]);

    DistanceToFacet distance_to_facet[24];

    Integer facet_index = -1;
    ENUMERATE_FACE(iface, cell.faces())
    {
      Face face = *iface;

      for (Integer i = 0; i < 4; i++)
      {
        facet_index++;

        Integer first_pos_node = (scan_order[iface.index()] ? ((i == 3) ? 0 : i+1) : i);
        Integer second_pos_node = (scan_order[iface.index()] ? i : ((i == 3) ? 0 : i+1));

        Node first_node = face.node(first_pos_node);
        Node second_node = face.node(second_pos_node);

        MC_Vector point0(m_coord_cm[first_node]);
        MC_Vector point1(m_coord_cm[second_node]);
        MC_Vector point2(m_coord_mid_cm[iface]);

        distance_to_facet[facet_index].distance = PhysicalConstants::_hugeDouble;

        MC_General_Plane plane(point0, point1, point2);

        Real facet_normal_dot_direction_cosine =
            (plane.A * m_particle_dir_cos[particle][MD_DirA] +
            plane.B * m_particle_dir_cos[particle][MD_DirB] +
            plane.C * m_particle_dir_cos[particle][MD_DirG]);

        // Consider only those facets whose outer normals have
        // a positive dot product with the direction cosine.
        // I.e. the particle is LEAVING the cell.
        if (facet_normal_dot_direction_cosine <= 0.0) continue;

        Real t = distanceToSegmentFacet(
            plane_tolerance,
            facet_normal_dot_direction_cosine, plane.A, plane.B, plane.C, plane.D,
            point0, point1, point2,
            particle, false);

        distance_to_facet[facet_index].distance = t;
      }
    }


    nearest_facet = findNearestFacet(
      particle,
      iteration, move_factor,
      distance_to_facet,
      retry);
  }

  if (nearest_facet.distance_to_facet < 0) 
  {
    nearest_facet.distance_to_facet = 0; 
  }

  if (nearest_facet.distance_to_facet >= PhysicalConstants::_hugeDouble)
  {
    qs_assert(false);
  }
  return nearest_facet;
}

/**
 * @brief 
 * 
 * @param plane_tolerance 
 * @param facet_normal_dot_direction_cosine 
 * @param A 
 * @param B 
 * @param C 
 * @param D 
 * @param facet_coords0 
 * @param facet_coords1 
 * @param facet_coords2 
 * @param particle 
 * @param allow_enter 
 * @return Real 
 */
Real TrackingMCModule::
distanceToSegmentFacet( Real plane_tolerance,
                        Real facet_normal_dot_direction_cosine,
                        Real A, Real B, Real C, Real D,
                        const MC_Vector &facet_coords0,
                        const MC_Vector &facet_coords1,
                        const MC_Vector &facet_coords2,
                        Particle particle,
                        bool allow_enter)
{
  Real boundingBox_tolerance = 1e-9;
  Real numerator = -1.0*(A * m_particle_coord[particle][MD_DirX] +
                            B * m_particle_coord[particle][MD_DirY] +
                            C * m_particle_coord[particle][MD_DirZ] +
                            D);

  /* Plane equation: numerator = -P(x,y,z) = -(Ax + By + Cz + D)
      if: numerator < -1e-8*length(x,y,z)   too negative!
      if: numerator < 0 && numerator^2 > ( 1e-8*length(x,y,z) )^2   too negative!
      reverse inequality since squaring function is decreasing for negative inputs.
      If numerator is just SLIGHTLY negative, then the particle is just outside of the face */

  // Filter out too negative distances
  if (!allow_enter && numerator < 0.0 && numerator * numerator > plane_tolerance)
  {
    return PhysicalConstants::_hugeDouble;
  }

  // we have to restrict the solution to within the triangular face
  Real distance = numerator / facet_normal_dot_direction_cosine;

  // see if the intersection point of the ray and the plane is within the triangular facet
  MC_Vector intersection_pt;
  intersection_pt.x = m_particle_coord[particle][MD_DirX] + distance * m_particle_dir_cos[particle][MD_DirA];
  intersection_pt.y = m_particle_coord[particle][MD_DirY] + distance * m_particle_dir_cos[particle][MD_DirB];
  intersection_pt.z = m_particle_coord[particle][MD_DirZ] + distance * m_particle_dir_cos[particle][MD_DirG];

  // if the point is completely below the triangle, it is not in the triangle
#define IF_POINT_BELOW_CONTINUE(axis)                                    \
  if ( facet_coords0.axis > intersection_pt.axis + boundingBox_tolerance&& \
        facet_coords1.axis > intersection_pt.axis + boundingBox_tolerance && \
        facet_coords2.axis > intersection_pt.axis + boundingBox_tolerance ) \
    { return PhysicalConstants::_hugeDouble; }

#define IF_POINT_ABOVE_CONTINUE(axis)                                    \
  if ( facet_coords0.axis < intersection_pt.axis - boundingBox_tolerance && \
        facet_coords1.axis < intersection_pt.axis - boundingBox_tolerance && \
        facet_coords2.axis < intersection_pt.axis - boundingBox_tolerance ) \
    { return PhysicalConstants::_hugeDouble; }

  // Is the intersection point inside the triangular facet?  Project to 2D and see.

  // A^2 + B^2 + C^2 = 1, so max(|A|,|B|,|C|) >= 1/sqrt(3) = 0.577
  // (all coefficients can't be small)
  Real cross0 = 0, cross1 = 0, cross2 = 0;
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

  Real cross_tol = 1e-9 * std::abs(cross0 + cross1 + cross2);  // cross product tolerance

  if ( (cross0 > -cross_tol && cross1 > -cross_tol && cross2 > -cross_tol) ||
      (cross0 <  cross_tol && cross1 <  cross_tol && cross2 <  cross_tol) )
  {
    return distance;
  }
  return PhysicalConstants::_hugeDouble;
}

/**
 * @brief Méthode permettant de trouver la facet la plus proche de la particule.
 * Cette méthode essaye de trouver la facet la plus proche en appelant nearestFacet.
 * Si on ne trouve pas, la méthode déplace la particule et rééssaye.
 * 
 */
NearestFacet TrackingMCModule::
findNearestFacet(Particle particle,
                  Integer &iteration, // input/output
                  Real &move_factor, // input/output
                  DistanceToFacet *distance_to_facet,
                  Integer &retry /* output */ )
{
  NearestFacet nearest_facet = nearestFacet(distance_to_facet);

  const Integer max_allowed_segments = 10000000;

  retry = 0;

  if ( (nearest_facet.distance_to_facet == PhysicalConstants::_hugeDouble && move_factor > 0) ||
      ( m_particle_num_seg[particle] > max_allowed_segments && nearest_facet.distance_to_facet <= 0.0 ) )
  {
    error() << "Attention, peut-être problème de facet.";
    error() << (nearest_facet.distance_to_facet == PhysicalConstants::_hugeDouble) << " && " << (move_factor > 0)
            << " || " << (m_particle_num_seg[particle] > max_allowed_segments) << " && " << (nearest_facet.distance_to_facet <= 0.0);

    // Could not find a solution, so move the particle towards the center of the cell
    // and try again.
    Cell cell = particle.cell();

    m_particle_coord[particle][MD_DirX] += move_factor * ( m_coord_center[cell][MD_DirX] - m_particle_coord[particle][MD_DirX] );
    m_particle_coord[particle][MD_DirY] += move_factor * ( m_coord_center[cell][MD_DirY] - m_particle_coord[particle][MD_DirY] );
    m_particle_coord[particle][MD_DirZ] += move_factor * ( m_coord_center[cell][MD_DirZ] - m_particle_coord[particle][MD_DirZ] );

    iteration++;
    move_factor *= 2.0;

    if ( move_factor > 1.0e-2 )
        move_factor = 1.0e-2;

    Integer max_iterations = 10000;

    if ( iteration == max_iterations )
    {
      qs_assert(false); // If we start hitting this assertion we can
      // come up with a better mitigation plan. - dfr
      retry = 0;
    }
    else
      retry = 1;

    // Allow the distance to the current facet
    //location->facet = -1;

  }
  return nearest_facet;
}

/**
 * @brief Méthode permettant de trouver la facet la plus proche de la particule.
 */
NearestFacet TrackingMCModule::
nearestFacet( DistanceToFacet *distance_to_facet)
{
  NearestFacet nearest_facet;

  // largest negative distance (smallest magnitude, but negative)
  NearestFacet nearest_negative_facet;
  nearest_negative_facet.distance_to_facet = -PhysicalConstants::_hugeDouble;

  // Determine the facet that is closest to the specified coordinates.
  for (Integer facet_index = 0; facet_index < 24; facet_index++)
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

/**
 * @brief Méthode permettant de trouver le plus petit élément d'un UniqueArray.
 * 
 * @tparam T Le type de l'UniqueArray.
 * @param array L'UniqueArray.
 * @return Integer La position du plus petit élément.
 */
template<typename T>
Integer TrackingMCModule::
findMin(UniqueArray<T> array)
{
  Real min = array[0];
  Integer min_index = 0;

  for (Integer element_index = 1; element_index < array.size(); ++element_index)
  {
    if ( array[element_index] < min )
    {
      min = array[element_index];
      min_index = element_index;
    }
  }

  return min_index;
}

/**
 * @brief Méthode permettant de faire une rotation à une particule.
 * 
 * @param particle La particule à tourner.
 * @param sin_Theta 
 * @param cos_Theta 
 * @param sin_Phi 
 * @param cos_Phi 
 */
void TrackingMCModule::
rotate3DVector(Particle particle, Real sin_Theta, Real cos_Theta, Real sin_Phi, Real cos_Phi)
{
  // Calculate additional variables in the rotation matrix.
  Real cos_theta = m_particle_dir_cos[particle][MD_DirG];
  Real sin_theta = sqrt((1.0 - (cos_theta*cos_theta)));

  Real cos_phi;
  Real sin_phi;
  if (sin_theta < 1e-6) // Order of sqrt(PhysicalConstants::tiny_double)
  {
    cos_phi = 1.0; // assume phi  = 0.0;
    sin_phi = 0.0;
  }
  else
  {
    cos_phi = m_particle_dir_cos[particle][MD_DirA]/sin_theta;
    sin_phi = m_particle_dir_cos[particle][MD_DirB]/sin_theta;
  }

  // Calculate the rotated direction cosine
  m_particle_dir_cos[particle][MD_DirA] =  cos_theta*cos_phi*(sin_Theta*cos_Phi) - sin_phi*(sin_Theta*sin_Phi) + sin_theta*cos_phi*cos_Theta;
  m_particle_dir_cos[particle][MD_DirB] =  cos_theta*sin_phi*(sin_Theta*cos_Phi) + cos_phi*(sin_Theta*sin_Phi) + sin_theta*sin_phi*cos_Theta;
  m_particle_dir_cos[particle][MD_DirG] = -sin_theta        *(sin_Theta*cos_Phi) +                               cos_theta        *cos_Theta;
}
