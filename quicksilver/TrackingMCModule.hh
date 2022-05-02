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
#include <arccore/concurrency/Mutex.h>
#include "structEnum.hh"

#include "TrackingMC_axl.h"

#include "NuclearData.hh"
#include "MC_Vector.hh"

using namespace Arcane;
using namespace Arcane::Materials;


/**
 * @brief Module TrackingMC.
 * Module permettant de suivre les particules contenues dans la famille
 * jusqu'a leur temps de census (ou leur sortie de maillage).
 */
class TrackingMCModule : public ArcaneTrackingMCObject {
public:
  explicit TrackingMCModule(const ModuleBuildInfo &mbi) 
  : ArcaneTrackingMCObject(mbi)
  , m_particle_family(nullptr)
  , m_cartesian_mesh(nullptr)
  , m_material_mng(nullptr)
  {}

public:
  void initModule() override;
  void cycleTracking() override;
  void endModule() override;

  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

protected:
  IItemFamily* m_particle_family;

  ICartesianMesh* m_cartesian_mesh;
  Arcane::Materials::IMeshMaterialMng* m_material_mng;
  NuclearData* m_nuclearData;

  Int32UniqueArray m_local_ids_exit;
  //Int32UniqueArray m_local_ids_processed;

  Int32UniqueArray m_local_ids_extra;

  Int32UniqueArray m_local_ids_extra_srcP;   //Particle localId source
  Int32UniqueArray m_local_ids_extra_cellId; //Cell localId dst
  Int64UniqueArray m_local_ids_extra_gId;    //Futur globalId
  Int64UniqueArray m_local_ids_extra_rns;    //Futur RNS
  RealUniqueArray  m_local_ids_extra_energyOut;    //updateTrajectory
  RealUniqueArray  m_local_ids_extra_angleOut;    //updateTrajectory
  RealUniqueArray  m_local_ids_extra_energyOut_pSrc;    //updateTrajectory
  RealUniqueArray  m_local_ids_extra_angleOut_pSrc;    //updateTrajectory

  Int32UniqueArray m_local_ids_out;
  Int32UniqueArray m_rank_out;

  Int32UniqueArray m_local_ids_in;

  std::atomic<Int64> m_num_segments_a{0};
  std::atomic<Int64> m_escape_a{0};
  std::atomic<Int64> m_census_a{0};
  std::atomic<Int64> m_collision_a{0};
  std::atomic<Int64> m_scatter_a{0};
  std::atomic<Int64> m_fission_a{0};
  std::atomic<Int64> m_absorb_a{0};
  std::atomic<Int64> m_produce_a{0};
  std::atomic<Int64> m_end_a{0};

  GlobalMutex m_mutex_exit;
  GlobalMutex m_mutex_extra;
  GlobalMutex m_mutex_out;
  GlobalMutex m_mutex_flux;
  
protected:
  void tracking();
  void updateTallies();
  void initNuclearData();
  bool isInGeometry(Integer pos, Cell cell);
  void cycleTrackingGuts( Particle particle );
  void cycleTrackingFunction(Particle particle);
  void collisionEventSuite();
  void computeNextEvent(Particle particle);
  Integer collisionEvent(Particle particle);
  void facetCrossingEvent(Particle particle);
  void reflectParticle(Particle particle);
  void cloneParticles(Int32UniqueArray idsSrc, Int32UniqueArray idsNew, Int64UniqueArray rnsNew);
  void cloneParticle(Particle pSrc, Particle pNew, Int64 rns);
  void updateTrajectory( Real energy, Real angle, Particle particle );
  void computeCrossSection();
  void weightedMacroscopicCrossSection(Cell cell, Integer energyGroup);
  Real macroscopicCrossSection(Integer reactionIndex, Cell cell, Integer isoIndex, Integer energyGroup);
  NearestFacet getNearestFacet( Particle particle);
  Real distanceToSegmentFacet(Real plane_tolerance,
                              Real facet_normal_dot_direction_cosine,
                              Real A, Real B, Real C, Real D,
                              const MC_Vector &facet_coords0,
                              const MC_Vector &facet_coords1,
                              const MC_Vector &facet_coords2,
                              Particle particle,
                              bool allow_enter);
  NearestFacet findNearestFacet( Particle particle,
                                  Integer &iteration,
                                  Real &move_factor,
                                  DistanceToFacet *distance_to_facet,
                                  Integer &retry);
  NearestFacet nearestFacet( DistanceToFacet *distance_to_facet);
  template<typename T>
  Integer findMin(UniqueArray<T> array);
  void rotate3DVector(Particle particle, Real sin_Theta, Real cos_Theta, Real sin_Phi, Real cos_Phi);

};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_TRACKINGMC(TrackingMCModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
