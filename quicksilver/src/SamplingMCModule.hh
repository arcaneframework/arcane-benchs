// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* SamplingMCModule.hh                                         (C) 2000-2022 */
/*                                                                           */
/* Module de sampling QAMA                                                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "structEnum.hh"
#include <arcane/IItemFamily.h>
#include <arcane/IMesh.h>
#include <arcane/IParallelMng.h>
#include <arcane/IParticleFamily.h>
#include <arcane/ITimeLoopMng.h>
#include <arcane/ModuleBuildInfo.h>
#include <arcane/materials/MeshMaterialVariableRef.h>
#include <arccore/concurrency/Mutex.h>
#include "ISimpleTableOutput.hh"
#include <arcane/ServiceBuilder.h>

#include "SamplingMC_axl.h"

using namespace Arcane;

/**
 * @brief Module SamplingMC.
 * Module permettant de créer des particules et de controler la population.
 */
class SamplingMCModule : public ArcaneSamplingMCObject
{
 public:
  explicit SamplingMCModule(const ModuleBuildInfo& mbi)
  : ArcaneSamplingMCObject(mbi)
  , m_particle_family(nullptr)
  , m_timer(nullptr)
  , m_rr(0)
  , m_split(0)
  , m_start(0)
  {}

 public:
  void initModule() override;
  void cycleSampling() override;
  void cycleFinalize() override;
  void endModule() override;

  VersionInfo versionInfo() const override { return VersionInfo(1, 1, 0); }

 protected:
  IItemFamily* m_particle_family;
  ParticleVectorView m_processingView;
  Real m_source_particle_weight;

  Int64 m_start;
  std::atomic<Int64> m_source_a{ 0 };
  Int64 m_rr;
  Int64 m_split;

  GlobalMutex m_mutex;

  Timer* m_timer;

 protected:
  void clearCrossSectionCache();
  void setStatus();
  void sourceParticles();
  void populationControl();
  void initParticle(ParticleEnumerator p, const Int64& rns);
  void cloneParticles(Int32UniqueArray idsSrc, Int32UniqueArray idsNew, Int64UniqueArray rnsNew);
  void cloneParticle(Particle pSrc, Particle pNew, const Int64& rns);
  Real computeTetVolume(const Real3& v0_, const Real3& v1_, const Real3& v2_, const Real3& v3);
  void rouletteLowWeightParticles();
  void generate3DCoordinate(Particle p, VariableNodeReal3& node_coord);
  void sampleIsotropic(Particle p);
  Real getSpeedFromEnergy(Particle p);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_SAMPLINGMC(SamplingMCModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
