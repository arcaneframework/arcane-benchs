// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "InitMCModule.hh"
#include "NVTX_Range.hh"
#include "MC_RNG_State.hh"
#include "PhysicalConstants.hh"


#define MAX_PRODUCTION_SIZE 4

bool ordre_qs2[]  = {true, true, false, false, false, true};

Integer QS2ArcaneFacet2[] = {16, 17, 18, 19, 
                        7 , 6 , 5 , 4 , 
                        20, 23, 22, 21, 
                        8 , 9 , 10, 11, 
                        12, 13, 14, 15, 
                        3 , 2 , 1 , 0 };

Integer QS2ArcaneFace2[] = {4, 1, 5, 2, 3, 0};

Integer QS2ArcaneNode2[] = {0, 1, 2, 3,
                       0, 3, 2, 1,
                       1, 0, 3, 2,
                       0, 1, 2, 3,
                       0, 1, 2, 3,
                       0, 3, 2, 1};


void InitMCModule::
initModule()
{
  m_cartesian_mesh = ICartesianMesh::getReference(mesh(), true);
  m_material_mng = IMeshMaterialMng::getReference(mesh());

  m_particle_family = mesh()->findItemFamily("ArcaneParticles");
}

void InitMCModule::
endModule()
{
  ENUMERATE_CELL(icell, ownCells())
  {
    if(m_sourceRate[icell] != 10000000000) ARCANE_FATAL("Bizarre");
  }
}

void InitMCModule::
cycleInit()
{
  ENUMERATE_CELL(icell, ownCells())
  {
    if(m_sourceRate[icell] != 10000000000) ARCANE_FATAL("Bizarre");
  }
  clearCrossSectionCache();

  m_processingView = m_particle_family->view();

  sourceParticles();

  // Réduction ou augmentation du nombre de particules.
  populationControl(); // controls particle population

  // Roulette sur les particules avec faible poids.
  rouletteLowWeightParticles(); // Delete particles with low statistical weight
  updateTallies();
  ENUMERATE_CELL(icell, ownCells())
  {
    if(m_sourceRate[icell] != 10000000000) ARCANE_FATAL("Bizarre");
  }
}

void InitMCModule::
updateTallies()
{
  m_source = m_source_a;
  m_rr = m_rr_a;
  m_split = m_split_a;
  
  m_source_a = 0;
  m_rr_a = 0;
  m_split_a = 0;
}

void InitMCModule::
clearCrossSectionCache()
{
  ENUMERATE_CELL(icell, ownCells())
  {
    for(Integer i = 0; i < m_nGroups(); i++)
    {
      m_total[icell][i] = 0.0;
    }
  }
}


void InitMCModule::
sourceParticles()
{
  NVTX_Range range("MC_Source_Now");

  Real local_weight_particles = 0;

  ENUMERATE_CELL(icell, ownCells())
  {
    if(m_sourceRate[icell] != 10000000000) ARCANE_FATAL("Bizarre");
    Real cell_weight_particles = m_volume[icell] * m_sourceRate[icell] * m_global_deltat();
    local_weight_particles += cell_weight_particles;
  }

  Real total_weight_particles = 0;

  total_weight_particles = mesh()->parallelMng()->reduce(Parallel::ReduceSum, local_weight_particles);

  Int64 num_particles = options()->getNParticles();
  Real source_fraction = 0.1;
  Real source_particle_weight = total_weight_particles/(source_fraction * num_particles);
  // Store the source particle weight for later use.
  m_source_particle_weight = source_particle_weight;

  Int64 task_index = 0;
  Int64 particle_count = 0;


  ENUMERATE_CELL(icell, ownCells())
  {
    Real cell_weight_particles = m_volume[icell] * m_sourceRate[icell] * m_global_deltat();
    Real cell_num_particles_float = cell_weight_particles / source_particle_weight;
    particle_count += (Int64)cell_num_particles_float;
  }

  Int64UniqueArray uids(particle_count);
  Int32UniqueArray local_id_cells(particle_count);
  Int32UniqueArray particles_lid(particle_count);
  Int64UniqueArray rng(particle_count);
  Integer particle_index_g = 0;

  ENUMERATE_CELL(icell, ownCells())
  {
    Real cell_weight_particles = m_volume[icell] * m_sourceRate[icell] * m_global_deltat();
    Real cell_num_particles_float = cell_weight_particles / source_particle_weight;
    Integer cell_num_particles = (Integer)cell_num_particles_float;

    for ( Integer particle_index = 0; particle_index < cell_num_particles; particle_index++ )
    {
      Int64 random_number_seed;
      Int64 rns;
      Int64 id;

      random_number_seed = m_sourceTally[icell]; // TODO : Atomic
      m_sourceTally[icell]++; // TODO : Atomic

      //if(particle_index != random_number_seed) ARCANE_FATAL("aaaa");

      random_number_seed += (*icell).uniqueId().asInt64() * INT64_C(0x0100000000);

      rns = rngSpawn_Random_Number_Seed(&random_number_seed);
      id = random_number_seed;

      rng[particle_index_g] = rns;
      uids[particle_index_g] = id;
      local_id_cells[particle_index_g] = icell.localId();

      particle_index_g++;
    }
  }

  m_particle_family->toParticleFamily()->addParticles(uids, local_id_cells, particles_lid);
  m_particle_family->endUpdate();

  ParticleVectorView viewSrcP = m_particle_family->view(particles_lid);

  ENUMERATE_PARTICLE(ipartic, viewSrcP)
  {
    Particle p = (*ipartic);
    initParticle(p, rng[ipartic.index()]);

    generate3DCoordinate(p);
    sampleIsotropic(p);
    m_particleKinEne[p] = (m_eMax() - m_eMin())*
                            rngSample(&m_particleRNS[p]) + m_eMin();

    Real speed = getSpeedFromEnergy(p);

    m_particleVelocity[p][MD_DirX] = speed * m_particleDirCos[p][MD_DirA];
    m_particleVelocity[p][MD_DirY] = speed * m_particleDirCos[p][MD_DirB];
    m_particleVelocity[p][MD_DirZ] = speed * m_particleDirCos[p][MD_DirG];

    m_particleTask[p] = task_index;
    m_particleWeight[p] = source_particle_weight;

    Real randomNumber = rngSample(&m_particleRNS[p]);
    m_particleNumMeanFreeP[p] = -1.0*std::log(randomNumber);

    randomNumber = rngSample(&m_particleRNS[p]);
    m_particleTimeCensus[p] = m_global_deltat() * randomNumber;

    m_source_a++;
  }

  m_processingView = m_particle_family->view();
}


void InitMCModule::
populationControl()
{
  NVTX_Range range("populationControl");

  Int64 targetNumParticles = options()->getNParticles();
  Int64 globalNumParticles = 0;
  Integer localNumParticles = m_processingView.size();
  
//     if (loadBalance)
//     {
//       // If we are parallel, we will have one domain per mpi processs.  The targetNumParticles is across
//       // all MPI processes, so we need to divide by the number or ranks to get the per-mpi-process number targetNumParticles
//       targetNumParticles = ceil((Real)targetNumParticles / mesh()->parallelMng()->commSize() );

//       //NO LONGER SPLITING VAULTS BY THREADS
// //        // If we are threaded, targetNumParticles should be divided by the number of threads (tasks) to balance
// //        // the particles across the thread level vaults.
// //        targetNumParticles = ceil((Real)targetNumParticles / (Real)monteCarlo->processor_info->num_tasks);
//     }
//     else
//     {
    globalNumParticles = mesh()->parallelMng()->reduce(Parallel::ReduceSum, localNumParticles);
  // }
    
  Real splitRRFactor = 1.0;
  // if (loadBalance)
  // {
  //     Integer currentNumParticles = localNumParticles;
  //     if (currentNumParticles != 0)
  //         splitRRFactor = (Real)targetNumParticles / (Real)currentNumParticles;
  //     else
  //         splitRRFactor = 1.0;
  // }
  // else
  // {
      splitRRFactor = (Real)targetNumParticles / (Real)globalNumParticles;
  // }

  // On augmente ou diminue la population selon splitRRFactor (si > 1, on augmente en splittant ; si < 1, on diminue en killant ou en augmentant le poids (rand))
  if (splitRRFactor != 1.0)  // no need to split if population is already correct.
    populationControlGuts(splitRRFactor, localNumParticles);
}


void InitMCModule::
populationControlGuts(const Real splitRRFactor, Int64 currentNumParticles)
{
  Int32UniqueArray supprP;
  Int64UniqueArray addIdP;
  Int32UniqueArray addCellIdP;
  Int32UniqueArray addSrcP;

  // March backwards through the vault so killed particles doesn't mess up the indexing
  ENUMERATE_PARTICLE(iparticle, m_processingView)
  {
    Particle particle = (*iparticle);
    Real randomNumber = rngSample(&m_particleRNS[iparticle]);
    if (splitRRFactor < 1)
    {
      if (randomNumber > splitRRFactor)
      {
        // Kill
        supprP.add(particle.localId());

        m_rr_a++; 
      }
      else
      {
        m_particleWeight[iparticle] /= splitRRFactor;
      }
    }
    else if (splitRRFactor > 1)
    {
      // Split
      Integer splitFactor = (Integer)floor(splitRRFactor);
      if (randomNumber > (splitRRFactor - splitFactor)) { splitFactor--; }

      m_particleWeight[iparticle] /= splitRRFactor;

      for (Integer splitFactorIndex = 0; splitFactorIndex < splitFactor; splitFactorIndex++)
      {
        m_split_a++;

        Int64 rns = rngSpawn_Random_Number_Seed(&m_particleRNS[iparticle]);
        addIdP.add(rns);
        addCellIdP.add(particle.cell().localId());
        addSrcP.add(particle.localId());
      }
    }
  }
  if (splitRRFactor < 1)
  {
    m_particle_family->toParticleFamily()->removeParticles(supprP);
    m_particle_family->toParticleFamily()->endUpdate();
    m_processingView = m_particle_family->view();
  }
  else if (splitRRFactor > 1)
  {
    Int32UniqueArray particles_lid(addIdP.size());
    m_particle_family->toParticleFamily()->addParticles(addIdP, addCellIdP, particles_lid);
    m_particle_family->toParticleFamily()->endUpdate();
    m_processingView = m_particle_family->view();

    copyParticles(addSrcP, particles_lid);
  }
}


void InitMCModule::
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

void InitMCModule::
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


void InitMCModule::
rouletteLowWeightParticles()
{
  NVTX_Range range("rouletteLowWeightParticles");

  const Real lowWeightCutoff = options()->getLowWeightCutoff();

  if (lowWeightCutoff > 0.0)
  {
    Int32UniqueArray supprP;

    // March backwards through the vault so killed particles don't mess up the indexing
    const Real source_particle_weight = m_source_particle_weight;
    const Real weightCutoff = lowWeightCutoff * source_particle_weight;

    ENUMERATE_PARTICLE(iparticle, m_processingView)
    {
      if (m_particleWeight[iparticle] <= weightCutoff)
      {
        Real randomNumber = rngSample(&m_particleRNS[iparticle]);
        if (randomNumber <= lowWeightCutoff)
        {
          // The particle history continues with an increased weight.
          m_particleWeight[iparticle] /= lowWeightCutoff;
        }
        else
        {
          // Kill
          supprP.add(iparticle.localId());
          m_rr_a++;
        } 
      }
    }
    m_particle_family->toParticleFamily()->removeParticles(supprP);
    m_particle_family->toParticleFamily()->endUpdate();
    m_processingView = m_particle_family->view();

  }
}


void InitMCModule::
initParticle(Particle p, Int64 rns)
{
  m_particleRNS[p] = rns;
  m_particleCoord[p][MD_DirX] = 0.0;
  m_particleCoord[p][MD_DirY] = 0.0;
  m_particleCoord[p][MD_DirZ] = 0.0;

  m_particleVelocity[p][MD_DirX] = 0.0;
  m_particleVelocity[p][MD_DirY] = 0.0;
  m_particleVelocity[p][MD_DirZ] = 0.0;

  m_particleDirCos[p][MD_DirA] = 0.0;
  m_particleDirCos[p][MD_DirB] = 0.0;
  m_particleDirCos[p][MD_DirG] = 0.0;

  m_particleKinEne[p] = 0.0;
  m_particleWeight[p] = 0.0;
  m_particleTimeCensus[p] = 0.0;
  m_particleTotalCrossSection[p] = 0.0;
  m_particleAge[p] = 0.0;
  m_particleNumMeanFreeP[p] = 0.0;
  m_particleMeanFreeP[p] = 0.0;
  m_particleSegPathLength[p] = 0.0;
  m_particleLastEvent[p] = Tally_Event::Census1;
  m_particleNumColl[p] = 0;
  m_particleNumSeg[p] = 0.0;
  m_particleTask[p] = 0;
  m_particleSpecies[p] = 0;
  m_particleBreed[p] = 0;
  m_particleEneGrp[p] = 0;
  m_particleFace[p] = 0;
  m_particleFacet[p] = 0;
  m_particleNormalDot[p] = 0.0;
}


void InitMCModule::
generate3DCoordinate(Particle p)
{
  Cell cell = p.cell();
  Int64* random_number_seed = &m_particleRNS[p];

  // Determine the cell-center nodal point coordinates.
  MC_Vector center(m_coordCenter[cell]);

  Integer num_facets = 24;

  Real random_number = rngSample(random_number_seed);
  Real which_volume = random_number * 6.0 * m_volume[cell];

  // Find the tet to sample from.
  Real current_volume = 0.0;
  Integer facet_index = -1;
  Node first_node;
  Node second_node;
  Face face;

  for(Integer i = 0; i < 6; i++)
  {
    face = cell.face(QS2ArcaneFace2[i]);

    for (Integer j = 0; j < 4; j++)
    {
      facet_index++;

      Integer first_pos_node = QS2ArcaneNode2[i*4 + j];
      Integer second_pos_node = QS2ArcaneNode2[i*4 + ((j == 3) ? 0 : j+1)];

      first_node = face.node(first_pos_node);
      second_node = face.node(second_pos_node);

      MC_Vector point0(m_coordCm[first_node]);
      MC_Vector point1(m_coordCm[second_node]);
      MC_Vector point2(m_coordMidCm[face]);

      Real subvolume = computeTetVolume(point0, point1, point2, center);
      current_volume += subvolume;

      if(current_volume >= which_volume) { break; }
    }
    if(current_volume >= which_volume) { break; }
  }

  // Sample from the tet.
  Real r1 = rngSample(random_number_seed);
  Real r2 = rngSample(random_number_seed);
  Real r3 = rngSample(random_number_seed);

  // Cut and fold cube into prism.
  if (r1 + r2 > 1.0)
  {
      r1 = 1.0 - r1;
      r2 = 1.0 - r2;
  }
  // Cut and fold prism into tetrahedron.
  if (r2 + r3 > 1.0)
  {
      Real tmp = r3;
      r3 = 1.0 - r1 - r2;
      r2 = 1.0 - tmp;
  }
  else if (r1 + r2 + r3 > 1.0)
  {
      Real tmp = r3;
      r3 = r1 + r2 + r3 - 1.0;
      r1 = 1.0 - r2 - tmp;
  }

  // numbers 1-4 are the barycentric coordinates of the random point.
  Real r4 = 1.0 - r1 - r2 - r3;

  MC_Vector point0(m_coordCm[first_node]);
  MC_Vector point1(m_coordCm[second_node]);
  MC_Vector point2(m_coordMidCm[face]);

  m_particleCoord[p][MD_DirX] = ( r4 * center.x + r1 * point0.x + r2 * point1.x + r3 * point2.x );
  m_particleCoord[p][MD_DirY] = ( r4 * center.y + r1 * point0.y + r2 * point1.y + r3 * point2.y );
  m_particleCoord[p][MD_DirZ] = ( r4 * center.z + r1 * point0.z + r2 * point1.z + r3 * point2.z );
}


///  \return 6 times the volume of the tet.
///
///  subtract v3 from v0, v1 and v2.  Then take the triple product of v0, v1 and v2.
Real InitMCModule::
computeTetVolume(const MC_Vector &v0_, const MC_Vector &v1_, const MC_Vector &v2_, const MC_Vector &v3)
{
  MC_Vector v0(v0_), v1(v1_), v2(v2_);

  v0.x -= v3.x; v0.y -= v3.y; v0.z -= v3.z;
  v1.x -= v3.x; v1.y -= v3.y; v1.z -= v3.z;
  v2.x -= v3.x; v2.y -= v3.y; v2.z -= v3.z;

  return
    v0.z*(v1.x*v2.y - v1.y*v2.x) +
    v0.y*(v1.z*v2.x - v1.x*v2.z) +
    v0.x*(v1.y*v2.z - v1.z*v2.y);
}


void InitMCModule::
sampleIsotropic(Particle p)
{
  m_particleDirCos[p][MD_DirG] = 1.0 - 2.0*(uint64_t)rngSample(&m_particleRNS[p]);
  Real sine_gamma  = sqrt((1.0 - (m_particleDirCos[p][MD_DirG]*m_particleDirCos[p][MD_DirG])));
  Real phi         = PhysicalConstants::_pi*(2.0*(uint64_t)rngSample(&m_particleRNS[p]) - 1.0);

  m_particleDirCos[p][MD_DirA]  = sine_gamma * cos(phi);
  m_particleDirCos[p][MD_DirB]  = sine_gamma * sin(phi);
}


Real InitMCModule::
getSpeedFromEnergy(Particle p)
{
  Real energy = m_particleKinEne[p];
  static const Real rest_mass_energy = PhysicalConstants::_neutronRestMassEnergy;
  static const Real speed_of_light  = PhysicalConstants::_speedOfLight;


  return speed_of_light * sqrt(energy * (energy + 2.0*(rest_mass_energy)) /
                                ((energy + rest_mass_energy) * (energy + rest_mass_energy)));
}
