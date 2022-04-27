#include "TrackingMCModule.hh"
#include <set>
#include <map>
#include "PhysicalConstants.hh"
#include "MC_RNG_State.hh"
#include "MC_Facet_Geometry.hh"
#define MAX_PRODUCTION_SIZE 4

bool ordre_qs3[]  = {true, true, false, false, false, true};


void TrackingMCModule::
initModule()
{
  m_cartesian_mesh = ICartesianMesh::getReference(mesh(), true);
  m_material_mng = IMeshMaterialMng::getReference(mesh());

  m_particle_family = mesh()->findItemFamily("ArcaneParticles");
  m_particle_family->setHasUniqueIdMap(false);

  IParticleExchanger* pe = options()->particleExchanger();
  pe->initialize(m_particle_family);

  initNuclearData(); // Configuration des materiaux.
}

void TrackingMCModule::
endModule()
{

}


/// Initializes both the NuclearData and the MaterialDatabase.  These
/// two structures are inherently linked since the isotopeGids stored in
/// the MaterialDatabase must correspond to the isotope indices in the
/// NuclearData.
void TrackingMCModule::
initNuclearData()
{
  Integer nb_cross_section = options()->cross_section().size();

  std::map<String, Polynomial> crossSection;
  std::map<String, Real> crossSection2; // TODO : Tres Tres Moche


  for( Integer i = 0; i < nb_cross_section; i++ )
  {
    String cs_name = options()->cross_section[i].getName();

    Real aa = options()->cross_section[i].getA();
    Real bb = options()->cross_section[i].getB();
    Real cc = options()->cross_section[i].getC();
    Real dd = options()->cross_section[i].getD();
    Real ee = options()->cross_section[i].getE();
    Real nuBar = options()->cross_section[i].getNuBar();

    crossSection.insert(std::make_pair(cs_name, Polynomial(aa, bb, cc, dd, ee)));
    crossSection2.insert(std::make_pair(cs_name, nuBar));
  }

  // TODO : Si nb_cross_section == 0.

  Integer num_materials = options()->material().size();
  Integer num_isotopes = 0;

  for( Integer i = 0; i < num_materials; i++ )
  {
    num_isotopes += options()->material[i].getNIsotopes();
  }

  // TODO : Si num_materials == 0.

  m_nuclearData = new NuclearData(m_nGroups(), m_eMin(), m_eMax());
  
  m_nuclearData->_isotopes.reserve( num_isotopes );


  for( Integer i = 0; i < num_materials; i++ )
  {
    String mat_name = options()->material[i].getName();
    MeshMaterialInfo* newMaterial = m_material_mng->registerMaterialInfo(mat_name);
    info() << "Creation du materiau : " << mat_name;
    MeshEnvironmentBuildInfo ebi1(mat_name);
    ebi1.addMaterial(mat_name);
    m_material_mng->createEnvironment(ebi1);
  }
  m_material_mng->endCreate();

  Integer nb_geometry = options()->geometry().size();

  Int32UniqueArray localIdMat[nb_geometry];

  ENUMERATE_CELL(icell, ownCells())
  {
    for( Integer i = 0; i < nb_geometry; i++ )
    {
      if(isInGeometry(i, (*icell)))
      {
        localIdMat[i].add(icell.localId());
        break;
      }
    }
  }

  ConstArrayView<IMeshMaterial*> materials = m_material_mng->materials();
  info() << "Size materiau : " << materials.size();
  {
    MeshMaterialModifier modifier(m_material_mng);
    for( Integer i = 0; i < nb_geometry; i++ )
    {
      String materialName = options()->geometry[i].getMaterial();

      for (Integer j = 0; j < num_materials; j++)
      {
        info() << "Compare materiau : " << materials[j]->name() << " et " << materialName+"_"+materialName;

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
    String fissionCrossSection = options()->material[i].getFissionCrossSection();
    String scatteringCrossSection = options()->material[i].getScatteringCrossSection();
    String absorptionCrossSection = options()->material[i].getAbsorptionCrossSection();
    Real totalCrossSection = options()->material[i].getTotalCrossSection();
    Real fissionCrossSectionRatio = options()->material[i].getFissionCrossSectionRatio();
    Real scatteringCrossSectionRatio = options()->material[i].getScatteringCrossSectionRatio();
    Real absorptionCrossSectionRatio = options()->material[i].getAbsorptionCrossSectionRatio();

    Real nuBar = crossSection2.at(options()->material[i].getFissionCrossSection().localstr());

    ENUMERATE_MATCELL(icell, materials[i])
    {
      m_mass[icell] = mass;
      m_sourceRate[icell] = sourceRate;
    }
    m_isoGid.resize(nIsotopes);
    m_atomFraction.resize(nIsotopes);


    for (Integer iIso=0; iIso<nIsotopes; ++iIso)
    {
      Integer isotopeGid = m_nuclearData->addIsotope(
        nReactions,
        crossSection.at(fissionCrossSection),
        crossSection.at(scatteringCrossSection),
        crossSection.at(absorptionCrossSection),
        nuBar,
        totalCrossSection,
        fissionCrossSectionRatio,
        scatteringCrossSectionRatio,
        absorptionCrossSectionRatio);
      
      // atomFraction for each isotope is 1/nIsotopes.  Treats all
      // isotopes as equally prevalent.
      ENUMERATE_MATCELL(icell, materials[i])
      {
        m_isoGid[icell][iIso] = isotopeGid;
        m_atomFraction[icell][iIso] = 1.0/nIsotopes;
      }
    }
  }
}


void TrackingMCModule::
cycleTracking()
{
  tracking();
  updateTallies();
}


void TrackingMCModule::
tracking()
{
  bool done = false;
  ParticleVectorView inView;

  // Ne sert qu'a debug maintenant.
  m_local_ids_processed.clear();

  ParticleVectorView m_processingView = m_particle_family->view();

  // A partir d'ici, toutes les particles de m_particle_family sont à suivre.
  pinfo() << "sizeProcessing " << m_processingView.size();

  Integer particle_count = 0; // Initialize count of num_particles processed
  
  while (!done)
  {

    ENUMERATE_PARTICLE(iparticle, m_processingView)
    {
      Particle particle = (*iparticle);

      if(iparticle.index() % 50000 == 0)
      {
        pinfo() << "--------";
        pinfo() << iparticle.index() << "/" << m_processingView.size();
        pinfo() << "m_local_ids_processed : " << m_local_ids_processed.size() << " m_local_ids_extra : " << m_local_ids_extra.size();
        pinfo() << "--------";
      }
      cycleTrackingGuts(particle);
    }
    particle_count += m_processingView.size();

    
    
    if(mesh()->parallelMng()->commSize() > 1)
    {
      ENUMERATE_PARTICLE(iparticle, inView)
      {
        Particle particle = (*iparticle);
        cycleTrackingGuts(particle);
      }
      particle_count += inView.size();
    }

    pinfo() << "========";
    pinfo() << "m_local_ids_exit : " << m_local_ids_exit.size() << " m_local_ids_extra_cellId : " << m_local_ids_extra_cellId.size()<< " m_local_ids_extra : " << m_local_ids_extra.size();
    pinfo() << "m_local_ids_processed : " << m_local_ids_processed.size() << " m_local_ids_exit : " << m_local_ids_exit.size();
    pinfo() << "========";

    

    m_particle_family->toParticleFamily()->removeParticles(m_local_ids_exit);
    // endUpdate fait par collisionEventSuite;
    m_local_ids_exit.clear();

    collisionEventSuite();

    if(mesh()->parallelMng()->commSize() > 1)
    {
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

        pinfo() << "/////////";
        pinfo() << "Nb Particles out : " << m_local_ids_out.size();
        pinfo() << "Nb Particles in : " << m_local_ids_in.size();
        pinfo() << "/////////";

        m_rank_out.clear();
        m_local_ids_out.clear();
      } while(m_local_ids_in.size() == 0 && m_local_ids_extra.size() == 0 && !done);

      inView = m_particle_family->view(m_local_ids_in);
      m_local_ids_in.clear();
    }

    else if(m_local_ids_extra.size() == 0)
    {
      done = true;
    }

    m_processingView = m_particle_family->view(m_local_ids_extra);
    m_local_ids_extra.clear();
  }
  info() << "particle_count : " << particle_count;
}

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
  m_numSegments = m_numSegments_a;
  
  m_absorb_a = 0;
  m_census_a = 0;
  m_escape_a = 0;
  m_collision_a = 0;
  m_fission_a = 0;
  m_produce_a = 0;
  m_scatter_a = 0;
  m_numSegments_a = 0;
}


// Returns true if the specified coordinate in inside the specified
// geometry.  False otherwise
bool TrackingMCModule::
isInGeometry(Integer pos, Cell cell)
{
  bool inside = false;
  switch (options()->geometry[pos].getShape())
  {
    case eShape::BRICK:
      {
        if ( (m_coordCenter[cell][MD_DirX] >= options()->geometry[pos].getXMin() && m_coordCenter[cell][MD_DirX] <= options()->geometry[pos].getXMax()) &&
            (m_coordCenter[cell][MD_DirY] >= options()->geometry[pos].getYMin() && m_coordCenter[cell][MD_DirY] <= options()->geometry[pos].getYMax()) &&
            (m_coordCenter[cell][MD_DirZ] >= options()->geometry[pos].getZMin() && m_coordCenter[cell][MD_DirZ] <= options()->geometry[pos].getZMax()) )
          inside = true;
      }
      break;
    case eShape::SPHERE:
      {
        MC_Vector center(options()->geometry[pos].getXCenter(), options()->geometry[pos].getYCenter(), options()->geometry[pos].getZCenter());
        MC_Vector rr(m_coordCenter[cell]);
        if ( (rr-center).Length() <= options()->geometry[pos].getRadius())
          inside = true;
      }

      break;
    default:
      qs_assert(false);
  }
  return inside;
}


void TrackingMCModule::
cycleTrackingGuts( Particle particle )
{
  if ( m_particleTimeCensus[particle] <= 0.0 )
  {
      m_particleTimeCensus[particle] += m_global_deltat();
  }

  // Age
  if (m_particleAge[particle] < 0.0) 
  { 
    m_particleAge[particle] = 0.0;
  }

  //    Energy Group
  m_particleEneGrp[particle] = m_nuclearData->getEnergyGroup(m_particleKinEne[particle]);

  // set the particle.task to the index of the processed vault the particle will census into.
  m_particleTask[particle] = 0;//processed_vault;

  // loop over this particle until we cannot do anything more with it on this processor
  cycleTrackingFunction(particle);
  // if(particle_index < 11)
  // {
  // cout << "particle.identifier : " << mc_particle.identifier << endl;
  // cout << mc_particle.coordinate.x << " x " << mc_particle.coordinate.y << " x " << mc_particle.coordinate.z << endl;
  // }
  // else
  // {
  // exit(123);
  // }


  //Make sure this particle is marked as completed
  m_particleSpecies[particle] = -1;
}

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

      // Collision ou Census ou Facet crossing
      Segment_Outcome_type segment_outcome = computeNextEvent(particle);
      m_numSegments_a++;

      

      m_particleNumSeg[particle] += 1.;  /* Track the number of segments this particle has
                                          undergone this cycle on all processes. */
      switch (segment_outcome) {
      case Segment_Outcome_type::Collision:
        {
          // The particle undergoes a collision event producing:
          //   (0) Other-than-one same-species secondary particle, or
          //   (1) Exactly one same-species secondary particle.
          switch (collisionEvent(particle))
          {

          case 0:
            m_local_ids_exit.add(particle.localId());
            keepTrackingThisParticle = false;
            break;

          case 1:
            keepTrackingThisParticle = true;
            break;
          
          default:
            keepTrackingThisParticle = false;
            break;

          }
        }
        break;
  
      case Segment_Outcome_type::Facet_Crossing:
        {
          // The particle has reached a cell facet.
          Tally_Event facet_crossing_type = facetCrossingEvent(particle);

          // if(particle_index == 11 && DEBUG_compt < 10)
          // {
          //   info() << "   Passe ici !" << facet_crossing_type;
          // }
          if (facet_crossing_type == Tally_Event::Facet_Crossing_Transit_Exit)
          {
              keepTrackingThisParticle = true;  // Transit Event
          }
          else if (facet_crossing_type == Tally_Event::Facet_Crossing_Escape)
          {
              m_escape_a++;
              m_particleLastEvent[particle] = Tally_Event::Facet_Crossing_Escape;
              m_particleSpecies[particle] = -1;
              keepTrackingThisParticle = false;
              m_local_ids_exit.add(particle.localId());

          }
          else if (facet_crossing_type == Tally_Event::Facet_Crossing_Reflection)
          {
              reflectParticle(particle);

              keepTrackingThisParticle = true;
          }
          else
          {
              // Enters an adjacent cell in an off-processor domain.
              //mc_particle.species = -1;
              keepTrackingThisParticle = false;
              // Pas de m_local_ids_exit car la particle sera retirée de la famille par ExchangeParticles.
          }
        }
        break;
  
      case Segment_Outcome_type::Census:
        {
          // The particle has reached the end of the time step.
          m_local_ids_processed.add(particle.localId());
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


void TrackingMCModule::
collisionEventSuite()
{
  Int32UniqueArray particles_lid(m_local_ids_extra_gId.size());
  m_particle_family->toParticleFamily()->addParticles(m_local_ids_extra_gId, m_local_ids_extra_cellId, particles_lid);
  m_particle_family->toParticleFamily()->endUpdate();

  copyParticles(m_local_ids_extra_srcP, particles_lid);

  ParticleVectorView viewP = m_particle_family->view(particles_lid);

  for(Integer i = 0; i < particles_lid.size(); i++)
  {
    Particle particle = Particle(viewP[i].internal());
    //if(particle.uniqueId().asInt64() != m_local_ids_extra_gId[i]) ARCANE_FATAL("Erreur ordre particle");
    updateTrajectory(m_local_ids_extra_energyOut[i], m_local_ids_extra_angleOut[i], particle);
  }

  viewP = m_particle_family->view(m_local_ids_extra);

  for(Integer i = 0; i < m_local_ids_extra.size(); i++)
  {
    Particle particle = Particle(viewP[i].internal());
    updateTrajectory(m_local_ids_extra_energyOut_pSrc[i], m_local_ids_extra_angleOut_pSrc[i], particle);
    m_particleEneGrp[particle] = m_nuclearData->getEnergyGroup(m_particleKinEne[particle]);
  }

  //Integer debug_size = m_local_ids_extra.size();
  //m_local_ids_extra.copy(particles_lid);

  //m_local_ids_extra.resize(m_local_ids_extra.size() + particles_lid.size());
  for(Integer i = 0; i < particles_lid.size(); i++)
  {
    m_local_ids_extra.add(particles_lid[i]);
  }

  //if(m_local_ids_extra.size() != debug_size + particles_lid.size()) 
  //{
  //  error() << "m_local_ids_extra.size() : " << m_local_ids_extra.size() << " debug_size : " << debug_size << " particles_lid.size() : " << particles_lid.size();
  //  ARCANE_FATAL("TODO A modifier");
  //}

  m_local_ids_extra_gId.clear();
  m_local_ids_extra_cellId.clear();
  m_local_ids_extra_srcP.clear();
  m_local_ids_extra_energyOut.clear();
  m_local_ids_extra_angleOut.clear();
  m_local_ids_extra_energyOut_pSrc.clear();
  m_local_ids_extra_angleOut_pSrc.clear();
}

Segment_Outcome_type TrackingMCModule::
computeNextEvent(Particle particle)
{
  // initialize distances to large number
  Integer number_of_events = 3;
  RealUniqueArray distance(3);
  distance[0] = distance[1] = distance[2] = 1e80; // +inf

  // Calculate the particle speed
  MC_Vector velo(m_particleVelocity[particle]);
  Real particle_speed = velo.Length();

  // Force collision if a census event narrowly preempts a collision
  Integer force_collision = 0 ;
  if ( m_particleNumMeanFreeP[particle] < 0.0 )
  {
    force_collision = 1 ;

    if ( m_particleNumMeanFreeP[particle] > -900.0 )
    {
      printf(" computeNextEvent: m_particleNumMeanFreeP[particle] > -900.0 \n");
    }

    m_particleNumMeanFreeP[particle] = PhysicalConstants::_smallDouble;
  }

  // Randomly determine the distance to the next collision
  // based upon the composition of the current cell.
  Real macroscopic_total_cross_section = weightedMacroscopicCrossSection(particle.cell(), m_particleEneGrp[particle]);
    
  // Cache the cross section
  m_particleTotalCrossSection[particle] = macroscopic_total_cross_section;
  if (macroscopic_total_cross_section == 0.0)
  {
    m_particleMeanFreeP[particle] = PhysicalConstants::_hugeDouble;
  }
  else
  {
    m_particleMeanFreeP[particle] = 1.0 / macroscopic_total_cross_section;
  }

  if ( m_particleNumMeanFreeP[particle] == 0.0)
  {
    // Sample the number of mean-free-paths remaining before
    // the next collision from an exponential distribution.
    Real random_number = rngSample(&m_particleRNS[particle]);

    m_particleNumMeanFreeP[particle] = -1.0*std::log(random_number);
  }

  // Calculate the distances to collision, nearest facet, and census.

  // Forced collisions do not need to move far.
  if (force_collision)
  {
    distance[Segment_Outcome_type::Collision] = PhysicalConstants::_smallDouble;
  }
  else
  {
    distance[Segment_Outcome_type::Collision] = m_particleNumMeanFreeP[particle] * m_particleMeanFreeP[particle];
  }

  // process census
  distance[Segment_Outcome_type::Census] = particle_speed*m_particleTimeCensus[particle];


  //  DEBUG  Turn off threshold for now
  Real distance_threshold = 10.0 * PhysicalConstants::_hugeDouble;
  // Get the current winning distance.
  Real current_best_distance = PhysicalConstants::_hugeDouble;


  bool new_segment =  (m_particleNumSeg[particle] == 0 ||
                        m_particleLastEvent[particle] == Tally_Event::Collision1);

  // Calculate the minimum distance to each facet of the cell.
  Nearest_Facet nearest_facet = getNearestFacet(particle, distance_threshold, current_best_distance, new_segment);

  m_particleNormalDot[particle] = nearest_facet.dot_product;

  distance[Segment_Outcome_type::Facet_Crossing] = nearest_facet.distance_to_facet;
  //info() << "Distance : " << nearest_facet.distance_to_facet << " " << nearest_facet.facet;
  //ARCANE_FATAL("aaa");

  // Get out of here if the tracker failed to bound this particle's volume.
  if (m_particleLastEvent[particle] == Tally_Event::Facet_Crossing_Tracking_Error)
  {
    return Segment_Outcome_type::Facet_Crossing;
  }

  // Calculate the minimum distance to the selected events.

  // Force a collision (if required).
  if ( force_collision == 1 )
  {
    distance[Segment_Outcome_type::Facet_Crossing] = PhysicalConstants::_hugeDouble;
    distance[Segment_Outcome_type::Census]         = PhysicalConstants::_hugeDouble;
    distance[Segment_Outcome_type::Collision]      = PhysicalConstants::_tinyDouble ;
  }

  // we choose our segment outcome here
  Segment_Outcome_type segment_outcome =
      (Segment_Outcome_type) findMin(distance);
  

  if (distance[segment_outcome] < 0)
  {
    // MC_Fatal_Jump( "Negative distances to events are NOT permitted!\n"
    //                 "identifier              = %" PRIu64 "\n"
    //                 "(Collision              = %g,\n"
    //                 " Facet Crossing         = %g,\n"
    //                 " Census                 = %g,\n",
    //                 mc_particle.identifier,
    //                 distance[Segment_Outcome_type::Collision],
    //                 distance[Segment_Outcome_type::Facet_Crossing],
    //                 distance[Segment_Outcome_type::Census]);
    qs_assert(false);
  }
  
  m_particleSegPathLength[particle] = distance[segment_outcome];

  m_particleNumMeanFreeP[particle] -= m_particleSegPathLength[particle] / m_particleMeanFreeP[particle];

  // Before using segment_outcome as an index, verify it is valid
  if (segment_outcome < 0 || segment_outcome >= Segment_Outcome_type::Max_Number)
  {
    // ( "segment_outcome '%d' is invalid\n", (Integer)segment_outcome );
    qs_assert(false);
  }

  Tally_Event SegmentOutcome_to_LastEvent[Segment_Outcome_type::Max_Number] =
  {
    Tally_Event::Collision1,
    Tally_Event::Facet_Crossing_Transit_Exit,
    Tally_Event::Census1,
  };

  m_particleLastEvent[particle] = SegmentOutcome_to_LastEvent[segment_outcome];

  // Set the segment path length to be the minimum of
  //   (i)   the distance to collision in the cell, or
  //   (ii)  the minimum distance to a facet of the cell, or
  //   (iii) the distance to census at the end of the time step
  if (segment_outcome == Segment_Outcome_type::Collision)
  {
    m_particleNumMeanFreeP[particle] = 0.0;
  }

  else if (segment_outcome == Segment_Outcome_type::Facet_Crossing)
  {
    m_particleFace[particle] = nearest_facet.facet / 4;
    m_particleFacet[particle] = nearest_facet.facet; // TODO : pos global ([0, 24[), voir pour mettre pos local ([0, 4[)
  }
  
  else if (segment_outcome == Segment_Outcome_type::Census)
  {
    m_particleTimeCensus[particle] = std::min(m_particleTimeCensus[particle], 0.0);
  }

  // If collision was forced, set mc_particle.num_mean_free_paths = 0
  // so that a new value is randomly selected on next pass.
  if (force_collision == 1)
  {
    m_particleNumMeanFreeP[particle] = 0.0;
  }

  // Do not perform any tallies if the segment path length is zero.
  //   This only introduces roundoff errors.
  if (m_particleSegPathLength[particle] == 0.0)
  {
    return segment_outcome;
  }

  // Move particle to end of segment, accounting for some physics processes along the segment.

  // Project the particle trajectory along the segment path length.

  //mc_particle.Move_Particle(mc_particle.direction_cosine, mc_particle.segment_path_length);
  m_particleCoord[particle][MD_DirX] += (m_particleDirCos[particle][MD_DirA] * m_particleSegPathLength[particle]);
  m_particleCoord[particle][MD_DirY] += (m_particleDirCos[particle][MD_DirB] * m_particleSegPathLength[particle]);
  m_particleCoord[particle][MD_DirZ] += (m_particleDirCos[particle][MD_DirG] * m_particleSegPathLength[particle]);

  Real segment_path_time = (m_particleSegPathLength[particle]/particle_speed);

  // Decrement the time to census and increment age.
  m_particleTimeCensus[particle] -= segment_path_time;
  m_particleAge[particle] += segment_path_time;

    //info() << "Drap" << mc_particle.age;
  // Ensure mc_particle.time_to_census is non-negative.
  if (m_particleTimeCensus[particle] < 0.0)
  {
    m_particleTimeCensus[particle] = 0.0;
  }

  // Accumulate the particle's contribution to the scalar flux.
  // TODO : Atomic
  m_scalarFluxTally[particle.cell()][m_particleEneGrp[particle]] += m_particleSegPathLength[particle] * m_particleWeight[particle];
  return segment_outcome;
}


Integer TrackingMCModule::
collisionEvent(Particle particle)
{
  Cell cell = particle.cell();

  //------------------------------------------------------------------------------------------------------------------
  //    Pick the isotope and reaction.
  //------------------------------------------------------------------------------------------------------------------
  Real randomNumber = rngSample(&m_particleRNS[particle]);
  Real totalCrossSection = m_particleTotalCrossSection[particle];
  Real currentCrossSection = totalCrossSection * randomNumber;
  Integer selectedIso = -1;
  Integer selectedUniqueNumber = -1;
  Integer selectedReact = -1;
  Integer numIsos = m_isoGid[cell].size();
  
  for (Integer isoIndex = 0; isoIndex < numIsos && currentCrossSection >= 0; isoIndex++)
  {
    Integer uniqueNumber = m_isoGid[cell][isoIndex];
    Integer numReacts = m_nuclearData->getNumberReactions(uniqueNumber);
    for (Integer reactIndex = 0; reactIndex < numReacts; reactIndex++)
    {
      currentCrossSection -= macroscopicCrossSection(reactIndex, cell,
                isoIndex, m_particleEneGrp[particle]);
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
  //Real mat_mass = monteCarlo->_materialDatabase->_mat[globalMatIndex]._mass;
  Real mat_mass = m_mass[cell];

  m_nuclearData->_isotopes[selectedUniqueNumber]._species[0]._reactions[selectedReact].sampleCollision(
    m_particleKinEne[particle], mat_mass, &energyOut[0], &angleOut[0], nOut, &(m_particleRNS[particle]), MAX_PRODUCTION_SIZE );

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
    m_particleEneGrp[particle] = m_nuclearData->getEnergyGroup(m_particleKinEne[particle]);
  }

  else
  {
    for (Integer secondaryIndex = 1; secondaryIndex < nOut; secondaryIndex++)
    {
      Int64 rns = rngSpawn_Random_Number_Seed(&m_particleRNS[particle]);
      m_local_ids_extra_gId.add(rns);
      m_local_ids_extra_cellId.add(particle.cell().localId());
      m_local_ids_extra_srcP.add(particle.localId());
      m_local_ids_extra_energyOut.add(energyOut[secondaryIndex]);
      m_local_ids_extra_angleOut.add(angleOut[secondaryIndex]);
    }

    // If a fission reaction produces secondary particles we also add the original
    // particle to the "extras" that we will handle later.  This avoids the 
    // possibility of a particle doing multiple fission reactions in a single
    // kernel invocation and overflowing the extra storage with secondary particles.
    m_local_ids_extra.add(particle.localId());
    m_local_ids_extra_energyOut_pSrc.add(energyOut[0]);
    m_local_ids_extra_angleOut_pSrc.add(angleOut[0]);
  }

  return nOut;
}


Tally_Event TrackingMCModule::
facetCrossingEvent(Particle particle)
{
  Face face = particle.cell().face(m_particleFace[particle]);

  if ( m_boundaryCond[face] == Face_Adjacency_Event::Transit_On_Processor )
  {
    // The particle will enter into an adjacent cell.
    Cell cell = face.frontCell();

    if(cell == particle.cell())
    {
      cell = face.backCell();
    }
    
    m_particleFacet[particle]     = (m_particleFacet[particle] < 12 ? m_particleFacet[particle] + 12 : m_particleFacet[particle] - 12);
    m_particleLastEvent[particle] = Tally_Event::Facet_Crossing_Transit_Exit;

    m_particle_family->toParticleFamily()->setParticleCell(particle, cell);
  }
  else if ( m_boundaryCond[face] == Face_Adjacency_Event::Boundary_Escape )
  {
    // The particle will escape across the system boundary.
    m_particleLastEvent[particle] = Tally_Event::Facet_Crossing_Escape;
  }
  else if ( m_boundaryCond[face] == Face_Adjacency_Event::Boundary_Reflection )
  {
    // The particle will reflect off of the system boundary.
    m_particleLastEvent[particle] = Tally_Event::Facet_Crossing_Reflection;
  }
  else if ( m_boundaryCond[face] == Face_Adjacency_Event::Transit_Off_Processor )
  {
    // The particle will enter into an adjacent cell on a spatial neighbor.
    
    Cell cell = face.frontCell();

    if(cell == particle.cell())
    {
      cell = face.backCell();
    }
    
    m_particle_family->toParticleFamily()->setParticleCell(particle, cell);

    m_particleFacet[particle]     = (m_particleFacet[particle] < 12 ? m_particleFacet[particle] + 12 : m_particleFacet[particle] - 12);
    m_particleLastEvent[particle] = Tally_Event::Facet_Crossing_Communication;

    m_local_ids_out.add(particle.localId());
    m_rank_out.add(cell.owner());
  }

  return (Tally_Event) m_particleLastEvent[particle];
}


void TrackingMCModule::
reflectParticle(Particle particle)
{
    Integer facet = m_particleFacet[particle] % 4;

    Cell cell = particle.cell();
    Face face = cell.face(m_particleFace[particle]);

    Integer first_pos_node = (ordre_qs3[m_indexArc[face]] ? ((facet == 3) ? 0 : facet+1) : facet);
    Integer second_pos_node = (ordre_qs3[m_indexArc[face]] ? facet : ((facet == 3) ? 0 : facet+1));

    Node first_node = face.node(first_pos_node);
    Node second_node = face.node(second_pos_node);

    MC_Vector point0(m_coordCm[first_node]);
    MC_Vector point1(m_coordCm[second_node]);
    MC_Vector point2(m_coordMidCm[face]);

    MC_General_Plane plane(point0, point1, point2);

    MC_Vector facet_normal(plane.A, plane.B, plane.C);


    Real dot = 2.0*( m_particleDirCos[particle][MD_DirA] * facet_normal.x +
                       m_particleDirCos[particle][MD_DirB] * facet_normal.y +
                       m_particleDirCos[particle][MD_DirG] * facet_normal.z );

    if ( dot > 0 ) // do not reflect a particle that is ALREADY pointing inward
    {
        // reflect the particle
        m_particleDirCos[particle][MD_DirA] -= dot * facet_normal.x;
        m_particleDirCos[particle][MD_DirB] -= dot * facet_normal.y;
        m_particleDirCos[particle][MD_DirG] -= dot * facet_normal.z;
    }

    // Calculate the reflected, velocity components.
    MC_Vector velo(m_particleVelocity[particle]);
    Real particle_speed = velo.Length();
    m_particleVelocity[particle][MD_DirX] = particle_speed * m_particleDirCos[particle][MD_DirA];
    m_particleVelocity[particle][MD_DirY] = particle_speed * m_particleDirCos[particle][MD_DirB];
    m_particleVelocity[particle][MD_DirZ] = particle_speed * m_particleDirCos[particle][MD_DirG];
}


void TrackingMCModule::
copyParticles(Int32UniqueArray idsSrc, Int32UniqueArray idsNew)
{
  ParticleVectorView viewSrcP = m_particle_family->view(idsSrc);
  ParticleVectorView viewNewP = m_particle_family->view(idsNew);
  for (Integer i = 0; i < idsSrc.size(); i++)
  {
    Particle pSrc(viewSrcP[i].internal());
    Particle pNew(viewNewP[i].internal());
    copyParticle(pSrc, pNew);
  }
}

void TrackingMCModule::
copyParticle(Particle pSrc, Particle pNew)
{
  m_particleRNS[pNew] = pNew.uniqueId().asInt64();

  m_particleCoord[pNew][MD_DirX] = m_particleCoord[pSrc][MD_DirX];
  m_particleCoord[pNew][MD_DirY] = m_particleCoord[pSrc][MD_DirY];
  m_particleCoord[pNew][MD_DirZ] = m_particleCoord[pSrc][MD_DirZ];

  m_particleVelocity[pNew][MD_DirX] = m_particleVelocity[pSrc][MD_DirX];
  m_particleVelocity[pNew][MD_DirY] = m_particleVelocity[pSrc][MD_DirY];
  m_particleVelocity[pNew][MD_DirZ] = m_particleVelocity[pSrc][MD_DirZ];

  m_particleDirCos[pNew][MD_DirA] = m_particleDirCos[pSrc][MD_DirA];
  m_particleDirCos[pNew][MD_DirB] = m_particleDirCos[pSrc][MD_DirB];
  m_particleDirCos[pNew][MD_DirG] = m_particleDirCos[pSrc][MD_DirG];

  m_particleKinEne[pNew] = m_particleKinEne[pSrc];
  m_particleWeight[pNew] = m_particleWeight[pSrc];
  m_particleTimeCensus[pNew] = m_particleTimeCensus[pSrc];
  m_particleTotalCrossSection[pNew] = m_particleTotalCrossSection[pSrc];
  m_particleAge[pNew] = m_particleAge[pSrc];
  m_particleNumMeanFreeP[pNew] = m_particleNumMeanFreeP[pSrc];
  m_particleMeanFreeP[pNew] = m_particleMeanFreeP[pSrc];
  m_particleSegPathLength[pNew] = m_particleSegPathLength[pSrc];
  m_particleLastEvent[pNew] = m_particleLastEvent[pSrc];
  m_particleNumColl[pNew] = m_particleNumColl[pSrc];
  m_particleNumSeg[pNew] = m_particleNumSeg[pSrc];
  m_particleTask[pNew] = m_particleTask[pSrc];
  m_particleSpecies[pNew] = m_particleSpecies[pSrc];
  m_particleBreed[pNew] = m_particleBreed[pSrc];
  m_particleEneGrp[pNew] = m_particleEneGrp[pSrc];
  m_particleFace[pNew] = m_particleFace[pSrc];
  m_particleFacet[pNew] = m_particleFacet[pSrc];
  m_particleNormalDot[pNew] = m_particleNormalDot[pSrc];
}

void TrackingMCModule::
updateTrajectory( Real energy, Real angle, Particle particle )
{
  m_particleKinEne[particle] = energy;
  Real cosTheta = angle;
  Real randomNumber = rngSample(&m_particleRNS[particle]);
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

  m_particleVelocity[particle][MD_DirX] = speed * m_particleDirCos[particle][MD_DirA];
  m_particleVelocity[particle][MD_DirY] = speed * m_particleDirCos[particle][MD_DirB];
  m_particleVelocity[particle][MD_DirZ] = speed * m_particleDirCos[particle][MD_DirG];

  randomNumber = rngSample(&m_particleRNS[particle]);
  m_particleNumMeanFreeP[particle] = -1.0*std::log(randomNumber);
}


Real TrackingMCModule::
weightedMacroscopicCrossSection(Cell cell, Integer energyGroup)
{
  Real precomputedCrossSection = m_total[cell][energyGroup];

  if (precomputedCrossSection > 0.0)
    return precomputedCrossSection;
  
  Integer nIsotopes = m_isoGid[cell].size();
  Real sum = 0.0;
  for (Integer isoIndex = 0; isoIndex < nIsotopes; isoIndex++)
  {
    sum += macroscopicCrossSection(-1, cell, isoIndex, energyGroup);
  }

  m_total[cell][energyGroup] = sum; // TODO Atomic

  return sum;
}

//----------------------------------------------------------------------------------------------------------------------
//  Routine MacroscopicCrossSection calculates the number-density-weighted macroscopic cross
//  section of a cell.
//
//  A reactionIndex of -1 means total cross section.
//----------------------------------------------------------------------------------------------------------------------
Real TrackingMCModule::
macroscopicCrossSection(Integer reactionIndex, Cell cell, Integer isoIndex, Integer energyGroup)
{
  // Initialize various data items.

  Real atomFraction = m_atomFraction[cell][isoIndex];

  Real microscopicCrossSection = 0.0;
  // The cell number density is the fraction of the atoms in cell
  // volume of this isotope.  We set this (elsewhere) to 1/nIsotopes.
  // This is a statement that we treat materials as if all of their
  // isotopes are present in equal amounts
  Real cellNumberDensity = m_cellNumberDensity[cell];

  Integer isotopeGid = m_isoGid[cell][isoIndex];
  if ( atomFraction == 0.0 || cellNumberDensity == 0.0) 
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

  return atomFraction * cellNumberDensity * microscopicCrossSection;
}

Nearest_Facet TrackingMCModule::
getNearestFacet( Particle particle,
                      Real distance_threshold,
                      Real current_best_distance,
                      bool new_segment)
{
  Nearest_Facet nearest_facet = computeFindNearestFacet(particle);

  if (nearest_facet.distance_to_facet < 0) {
    nearest_facet.distance_to_facet = 0; 
  }

  if (nearest_facet.distance_to_facet >= PhysicalConstants::_hugeDouble)
  {
    qs_assert(false);
  }

  return nearest_facet;
}

Nearest_Facet TrackingMCModule::
computeFindNearestFacet(Particle particle)
{
  Cell cell = particle.cell();
  MC_Vector *facet_coords[3];
  Integer iteration = 0;
  Real move_factor = 0.5 * PhysicalConstants::_smallDouble;

  while (true) // will break out when distance is found
  {
    // Determine the distance to each facet of the cell.
    // (1e-8 * Radius)^2
    Real plane_tolerance = 1e-16*(m_particleCoord[particle][MD_DirX]*m_particleCoord[particle][MD_DirX] +
                                    m_particleCoord[particle][MD_DirY]*m_particleCoord[particle][MD_DirY] +
                                    m_particleCoord[particle][MD_DirZ]*m_particleCoord[particle][MD_DirZ]);

    Distance_To_Facet distance_to_facet[24];

    Integer facet_index = -1;
    ENUMERATE_FACE(iface, cell.faces())
    {
      Face face = *iface;
      

      for (Integer i = 0; i < 4; i++)
      {
        facet_index++;
        //if((Integer)facet_index/4 != iface.index())
        //{
        //  ARCANE_FATAL("Erreur facet index");
        //}
        Integer first_pos_node = (ordre_qs3[iface.index()] ? ((i == 3) ? 0 : i+1) : i);
        Integer second_pos_node = (ordre_qs3[iface.index()] ? i : ((i == 3) ? 0 : i+1));

        Node first_node = face.node(first_pos_node);
        Node second_node = face.node(second_pos_node);

        MC_Vector point0(m_coordCm[first_node]);
        MC_Vector point1(m_coordCm[second_node]);
        MC_Vector point2(m_coordMidCm[iface]);

        distance_to_facet[facet_index].distance = PhysicalConstants::_hugeDouble;

        MC_General_Plane plane(point0, point1, point2);

        Real facet_normal_dot_direction_cosine =
            (plane.A * m_particleDirCos[particle][MD_DirA] +
            plane.B * m_particleDirCos[particle][MD_DirB] +
            plane.C * m_particleDirCos[particle][MD_DirG]);

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
    //ARCANE_FATAL("aaa");

    Integer retry = 0;

    Nearest_Facet nearest_facet = findNearestFacet(
      particle,
      iteration, move_factor,
      distance_to_facet,
      retry);


    if (!retry) return nearest_facet;
  }
}

Real TrackingMCModule::
distanceToSegmentFacet(Real plane_tolerance,
                                                     Real facet_normal_dot_direction_cosine,//=
                                                     Real A, Real B, Real C, Real D,//=
                                                     const MC_Vector &facet_coords0,//=
                                                     const MC_Vector &facet_coords1,//=
                                                     const MC_Vector &facet_coords2,//=
                                                     Particle particle,//=
                                                     bool allow_enter) //=
{
  Real boundingBox_tolerance = 1e-9;
  Real numerator = -1.0*(A * m_particleCoord[particle][MD_DirX] +
                            B * m_particleCoord[particle][MD_DirY] +
                            C * m_particleCoord[particle][MD_DirZ] +
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
  intersection_pt.x = m_particleCoord[particle][MD_DirX] + distance * m_particleDirCos[particle][MD_DirA];
  intersection_pt.y = m_particleCoord[particle][MD_DirY] + distance * m_particleDirCos[particle][MD_DirB];
  intersection_pt.z = m_particleCoord[particle][MD_DirZ] + distance * m_particleDirCos[particle][MD_DirG];

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

Nearest_Facet TrackingMCModule::
findNearestFacet(Particle particle,
                                  Integer &iteration, // input/output
                                  Real &move_factor, // input/output
                                  Distance_To_Facet *distance_to_facet,
                                  Integer &retry /* output */ )
{
  Nearest_Facet nearest_facet = nearestFacet(distance_to_facet);

  const Integer max_allowed_segments = 10000000;

  retry = 0;

  if ( (nearest_facet.distance_to_facet == PhysicalConstants::_hugeDouble && move_factor > 0) ||
      ( m_particleNumSeg[particle] > max_allowed_segments && nearest_facet.distance_to_facet <= 0.0 ) )
  {
    info() << "Attention, peut-être problème de facet.";
    // Could not find a solution, so move the particle towards the center of the cell
    // and try again.
    nearestFacet3DMoveParticle(particle, move_factor);
    iteration++;
    move_factor *= 2.0;

    if ( move_factor > 1.0e-2 )
        move_factor = 1.0e-2;

    Integer max_iterations = 10000;

    if ( iteration == max_iterations )
    {
      //info() << (nearest_facet.distance_to_facet == PhysicalConstants::_hugeDouble) << (move_factor > 0) << 
      //(mc_particle->num_segments > max_allowed_segments) << (nearest_facet.distance_to_facet <= 0.0);

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

Nearest_Facet TrackingMCModule::
nearestFacet( Distance_To_Facet *distance_to_facet)
{
  Nearest_Facet nearest_facet;

  // largest negative distance (smallest magnitude, but negative)
  Nearest_Facet nearest_negative_facet;
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

void TrackingMCModule::
nearestFacet3DMoveParticle( Particle particle, // input/output: move this coordinate
                            Real move_factor)  // input: multiplication factor for move
{
  Cell cell = particle.cell();

  m_particleCoord[particle][MD_DirX] += move_factor * ( m_coordCenter[cell][MD_DirX] - m_particleCoord[particle][MD_DirX] );
  m_particleCoord[particle][MD_DirY] += move_factor * ( m_coordCenter[cell][MD_DirY] - m_particleCoord[particle][MD_DirY] );
  m_particleCoord[particle][MD_DirZ] += move_factor * ( m_coordCenter[cell][MD_DirZ] - m_particleCoord[particle][MD_DirZ] );
}

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


void TrackingMCModule::
rotate3DVector(Particle particle, Real sin_Theta, Real cos_Theta, Real sin_Phi, Real cos_Phi)
{
  // Calculate additional variables in the rotation matrix.
  Real cos_theta = m_particleDirCos[particle][MD_DirG];
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
    cos_phi = m_particleDirCos[particle][MD_DirA]/sin_theta;
    sin_phi = m_particleDirCos[particle][MD_DirB]/sin_theta;
  }

  // Calculate the rotated direction cosine
  m_particleDirCos[particle][MD_DirA] =  cos_theta*cos_phi*(sin_Theta*cos_Phi) - sin_phi*(sin_Theta*sin_Phi) + sin_theta*cos_phi*cos_Theta;
  m_particleDirCos[particle][MD_DirB] =  cos_theta*sin_phi*(sin_Theta*cos_Phi) + cos_phi*(sin_Theta*sin_Phi) + sin_theta*sin_phi*cos_Theta;
  m_particleDirCos[particle][MD_DirG] = -sin_theta        *(sin_Theta*cos_Phi) +                               cos_theta        *cos_Theta;
}
