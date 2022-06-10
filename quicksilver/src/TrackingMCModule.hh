// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* TrackingMCModule.hh                                         (C) 2000-2022 */
/*                                                                           */
/* Module de tracking QAMA                                                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "structEnum.hh"
#include <arcane/IAsyncParticleExchanger.h>
#include <arcane/IItemFamily.h>
#include <arcane/IMesh.h>
#include <arcane/IParallelMng.h>
#include <arcane/IParticleExchanger.h>
#include <arcane/IParticleFamily.h>
#include <arcane/ITimeLoopMng.h>
#include <arcane/ModuleBuildInfo.h>
#include <arcane/Timer.h>
#include <arcane/materials/IMeshMaterialMng.h>
#include <arcane/materials/MaterialVariableBuildInfo.h>
#include <arcane/materials/MeshEnvironmentBuildInfo.h>
#include <arcane/materials/MeshMaterialModifier.h>
#include <arcane/materials/MeshMaterialVariableRef.h>
#include <arccore/concurrency/Mutex.h>
#include "ISimpleTableOutput.hh"
#include <arcane/ServiceBuilder.h>

#include "TrackingMC_axl.h"

#include "NuclearData.hh"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * @brief Structure simple permettant de stocker la distance d'une facet à une particule.
 */
struct DistanceToFacet
{
  Real distance;
  Integer facet;

  DistanceToFacet()
  : distance(0.0)
  , facet(0)
  {}
};

/**
 * @brief Module TrackingMC.
 * Module permettant de suivre les particules contenues dans la famille
 * jusqu'a leur temps de census (ou leur sortie de maillage).
 */
class TrackingMCModule : public ArcaneTrackingMCObject
{
 public:
  explicit TrackingMCModule(const ModuleBuildInfo& mbi)
  : ArcaneTrackingMCObject(mbi)
  , m_particle_family(nullptr)
  , m_material_mng(nullptr)
  , m_timer(nullptr)
  , m_nuclearData(nullptr)
  , m_exited_particles_local_ids(0)
  , m_extra_particles_local_ids(0)
  , m_extra_particles_particle_src(0)
  , m_extra_particles_cellid_dst(0)
  , m_extra_particles_global_id(0)
  , m_extra_particles_rns(0)
  , m_extra_particles_energy_out(0)
  , m_extra_particles_angle_out(0)
  , m_extra_particles_energy_out_particle_src(0)
  , m_extra_particles_angle_out_particle_src(0)
  , m_outgoing_particles_local_ids(0)
  , m_outgoing_particles_rank_to(0)
  , m_escape(0)
  , m_end(0)
  , m_incoming(0)
  , m_outgoing(0)
  {}

 public:
  void initModule() override;
  void cycleTracking() override;
  void cycleFinalize() override;
  void endModule() override;

  VersionInfo versionInfo() const override { return VersionInfo(1, 7, 0); }

 protected:
  IItemFamily* m_particle_family;

  IMeshMaterialMng* m_material_mng;
  NuclearData* m_nuclearData;
  Timer* m_timer;

  Int32UniqueArray m_exited_particles_local_ids;

  // Les particules "extra" sont les particules qui n'ont pas fini leur itération
  // (et donc necessite une autre sous-itération).
  Int32UniqueArray m_extra_particles_local_ids;
  Int32UniqueArray m_extra_particles_particle_src; //Particle localId source
  Int32UniqueArray m_extra_particles_cellid_dst; //Cell localId dst
  Int64UniqueArray m_extra_particles_global_id; //Futur globalId
  Int64UniqueArray m_extra_particles_rns; //Futur RNS
  RealUniqueArray m_extra_particles_energy_out; //updateTrajectory
  RealUniqueArray m_extra_particles_angle_out; //updateTrajectory
  RealUniqueArray m_extra_particles_energy_out_particle_src; //updateTrajectory
  RealUniqueArray m_extra_particles_angle_out_particle_src; //updateTrajectory

  Int32UniqueArray m_outgoing_particles_local_ids;
  Int32UniqueArray m_outgoing_particles_rank_to;

  std::atomic<Int64> m_num_segments_a{ 0 };
  Int64 m_escape;
  std::atomic<Int64> m_census_a{ 0 };
  std::atomic<Int64> m_collision_a{ 0 };
  std::atomic<Int64> m_scatter_a{ 0 };
  std::atomic<Int64> m_fission_a{ 0 };
  std::atomic<Int64> m_absorb_a{ 0 };
  std::atomic<Int64> m_produce_a{ 0 };
  Int64 m_end;

  Int64 m_incoming;
  Int64 m_outgoing;

  GlobalMutex m_mutex_exit;
  GlobalMutex m_mutex_extra;
  GlobalMutex m_mutex_out;
  GlobalMutex m_mutex_flux;
  GlobalMutex m_mutex_lb;

  bool m_do_loop_lb;

 protected:
  void tracking();
  void updateTallies();
  void initNuclearData();
  bool isInGeometry(const Integer& pos, Cell cell);
  void cycleTrackingGuts(Particle particle, VariableNodeReal3& node_coord);
  void cycleTrackingFunction(Particle particle, VariableNodeReal3& node_coord);
  void collisionEventSuite();
  void computeNextEvent(Particle particle, VariableNodeReal3& node_coord);
  Integer collisionEvent(Particle particle);
  void facetCrossingEvent(Particle particle);
  void reflectParticle(Particle particle, VariableNodeReal3& node_coord);
  void cloneParticles(Int32UniqueArray idsSrc, Int32UniqueArray idsNew, Int64UniqueArray rnsNew);
  void cloneParticle(Particle pSrc, Particle pNew, const Int64& rns);
  void updateTrajectory(const Real& energy, const Real& angle, Particle particle);
  void computeCrossSection();
  void weightedMacroscopicCrossSection(Cell cell, const Integer& energyGroup);
  Real macroscopicCrossSection(const Integer& reactionIndex,
                               const Real& cell_number_density,
                               const Real& atom_fraction,
                               const Integer& isotopeGid,
                               const Integer& isoIndex,
                               const Integer& energyGroup);
  DistanceToFacet getNearestFacet(Particle particle, VariableNodeReal3& node_coord);
  Real distanceToSegmentFacet(const Real& plane_tolerance,
                              const Real& facet_normal_dot_direction_cosine,
                              const Real& A, const Real& B, const Real& C, const Real& D,
                              const Real3& facet_coords0,
                              const Real3& facet_coords1,
                              const Real3& facet_coords2,
                              Particle particle,
                              bool allow_enter);
  DistanceToFacet findNearestFacet(Particle particle,
                                   Integer& iteration,
                                   Real& move_factor,
                                   DistanceToFacet* distance_to_facet,
                                   Integer& retry);
  DistanceToFacet nearestFacet(DistanceToFacet* distance_to_facet);
  template <typename T>
  Integer findMin(UniqueArray<T> array);
  void rotate3DVector(Particle particle, const Real& sin_Theta, const Real& cos_Theta, const Real& sin_Phi, const Real& cos_Phi);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_TRACKINGMC(TrackingMCModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
