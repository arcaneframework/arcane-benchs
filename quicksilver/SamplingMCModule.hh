// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include <arcane/IParticleExchanger.h>
#include <arcane/ITimeLoopMng.h>
#include <arcane/geometry/IGeometryMng.h>
#include <arcane/cartesianmesh/CellDirectionMng.h>
#include <arcane/cartesianmesh/ICartesianMesh.h>
#include <arcane/IParallelMng.h>
#include "arcane/ModuleBuildInfo.h"
#include "arcane/IMesh.h"
#include "arcane/IMeshModifier.h"
#include "arcane/IItemFamily.h"
#include "arcane/IParticleFamily.h"
#include "arcane/ItemVector.h"
#include "arcane/IParticleExchanger.h"
#include "arcane/IAsyncParticleExchanger.h"
#include "arcane/ItemPrinter.h"
#include "arcane/IExtraGhostParticlesBuilder.h"
#include "arcane/materials/IMeshMaterialMng.h"
#include <arcane/materials/IMeshMaterial.h>
#include "arcane/materials/IMeshEnvironment.h"
#include "arcane/materials/IMeshBlock.h"
#include "arcane/materials/MeshMaterialModifier.h"
#include "arcane/materials/MeshMaterialVariableRef.h"
#include "arcane/materials/MeshEnvironmentVariableRef.h"
#include "arcane/materials/MaterialVariableBuildInfo.h"
#include "arcane/materials/MeshEnvironmentBuildInfo.h"
#include "structEnum.hh"

#include "SamplingMC_axl.h"

#include "MC_Vector.hh"


using namespace Arcane;

/**
 * @brief Module SamplingMC.
 * Module permettant de créer des particules et de controler la population.
 */
class SamplingMCModule : 
public ArcaneSamplingMCObject {
public:
  explicit SamplingMCModule(const ModuleBuildInfo &mbi) : ArcaneSamplingMCObject(mbi) {}

public:
  void initModule() override;
  void cycleInit() override;
  void endModule() override;
  
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

protected:
  IItemFamily* m_particle_family;
  ParticleVectorView m_processingView;
  Real m_source_particle_weight;

  std::atomic<Int64> m_start_a{0};
  std::atomic<Int64> m_source_a{0};
  std::atomic<Int64> m_rr_a{0};
  std::atomic<Int64> m_split_a{0};


protected:
  void updateTallies();
  void clearCrossSectionCache();
  void sourceParticles();
  void populationControl();
  void cloneParticles(Int32UniqueArray idsSrc, Int32UniqueArray idsNew, Int64UniqueArray rnsNew);
  void cloneParticle(Particle pSrc, Particle pNew, Int64 rns);
  Real computeTetVolume(const MC_Vector &v0_, const MC_Vector &v1_, const MC_Vector &v2_, const MC_Vector &v3);
  void rouletteLowWeightParticles();
  void initParticle(Particle p, Int64 rns);
  void generate3DCoordinate(Particle p);
  void sampleIsotropic(Particle p);
  Real getSpeedFromEnergy(Particle p);

};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_SAMPLINGMC(SamplingMCModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
