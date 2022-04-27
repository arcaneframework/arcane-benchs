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

enum eShape{UNDEFINED, BRICK, SPHERE}; // TODO : A deplacer (doit être defini avant QS_axl.h !).
enum eBoundaryCondition{reflect, escape, octant}; // TODO : A deplacer (doit être defini avant QS_axl.h !).

#include "TrackingMC_axl.h"

#include "NuclearData.hh"
#include "MC_Vector.hh"

using namespace Arcane;
using namespace Arcane::Materials;

enum Face_Adjacency_Event
{
  Adjacency_Undefined = 0,
  Boundary_Escape,
  Boundary_Reflection,
  Transit_On_Processor,
  Transit_Off_Processor
};

enum CosDir
{
  MD_DirA = 0, // Alpha
  MD_DirB,     // Beta
  MD_DirG      // Gamma
};

enum Segment_Outcome_type
{
  Initialize                    = -1,
  Collision                     = 0,
  Facet_Crossing                = 1,
  Census                        = 2,
  Max_Number                    = 3
};
enum Tally_Event
{
  Collision1,
  Facet_Crossing_Transit_Exit,
  Census1,
  Facet_Crossing_Tracking_Error,
  Facet_Crossing_Escape,
  Facet_Crossing_Reflection,
  Facet_Crossing_Communication
};

struct Nearest_Facet
{
   Integer facet;
   Real distance_to_facet;
   Real dot_product;
   
   Nearest_Facet()
   : facet(0),
     distance_to_facet(1e80),
     dot_product(0.0)
   {}
};

struct Distance_To_Facet
{
    Real distance;
    Integer facet;
    Integer subfacet;
    
    Distance_To_Facet()
    : distance(0.0),
      facet(0),
      subfacet(0) 
    {}
};

/*!
 * \brief Module TrackingMC.
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

public:
  IItemFamily* m_particle_family;

  ICartesianMesh* m_cartesian_mesh;
  Arcane::Materials::IMeshMaterialMng* m_material_mng;

  Int32UniqueArray m_local_ids_exit;
  Int32UniqueArray m_local_ids_processed;

  Int32UniqueArray m_local_ids_extra;

  Int32UniqueArray m_local_ids_extra_srcP;   //Particle localId source
  Int32UniqueArray m_local_ids_extra_cellId; //Cell localId dst
  Int64UniqueArray m_local_ids_extra_gId;    //Futur globalId
  RealUniqueArray  m_local_ids_extra_energyOut;    //updateTrajectory
  RealUniqueArray  m_local_ids_extra_angleOut;    //updateTrajectory
  RealUniqueArray  m_local_ids_extra_energyOut_pSrc;    //updateTrajectory
  RealUniqueArray  m_local_ids_extra_angleOut_pSrc;    //updateTrajectory

  Int32UniqueArray m_local_ids_out;
  Int32UniqueArray m_rank_out;

  Int32UniqueArray m_local_ids_in;
  NuclearData* m_nuclearData;
  std::atomic<Int64> m_numSegments_a{0};
  std::atomic<Int64> m_escape_a{0};
  std::atomic<Int64> m_census_a{0};
  std::atomic<Int64> m_collision_a{0};
  std::atomic<Int64> m_scatter_a{0};
  std::atomic<Int64> m_fission_a{0};
  std::atomic<Int64> m_absorb_a{0};
  std::atomic<Int64> m_produce_a{0};
  
public:
  void tracking();
  void updateTallies();
  void initNuclearData();
  bool isInGeometry(Integer pos, Cell cell);
  void cycleTrackingGuts( Particle particle );
  void cycleTrackingFunction(Particle particle);
  void collisionEventSuite();
  Segment_Outcome_type computeNextEvent(Particle particle);
  Integer collisionEvent(Particle particle);
  Tally_Event facetCrossingEvent(Particle particle);
  void reflectParticle(Particle particle);
  void copyParticles(Int32UniqueArray idsSrc, Int32UniqueArray idsNew);
  void copyParticle(Particle pSrc, Particle pNew);
  void updateTrajectory( Real energy, Real angle, Particle particle );
  Real weightedMacroscopicCrossSection(Cell cell, Integer energyGroup);
  Real macroscopicCrossSection(Integer reactionIndex, Cell cell, Integer isoIndex, Integer energyGroup);
  Nearest_Facet getNearestFacet( Particle particle, Real distance_threshold,
                        Real current_best_distance,
                        bool new_segment);
  Nearest_Facet computeFindNearestFacet(Particle particle);
  Real distanceToSegmentFacet(Real plane_tolerance,
                                                      Real facet_normal_dot_direction_cosine,//=
                                                      Real A, Real B, Real C, Real D,//=
                                                      const MC_Vector &facet_coords0,//=
                                                      const MC_Vector &facet_coords1,//=
                                                      const MC_Vector &facet_coords2,//=
                                                      Particle particle,//=
                                                      bool allow_enter); //=
  Nearest_Facet findNearestFacet(Particle particle,
                                    Integer &iteration, // input/output
                                    Real &move_factor, // input/output
                                    Distance_To_Facet *distance_to_facet,
                                    Integer &retry /* output */ );
  Nearest_Facet nearestFacet( Distance_To_Facet *distance_to_facet);
  void nearestFacet3DMoveParticle( Particle particle, // input/output: move this coordinate
                              Real move_factor);
  template<typename T>
  Integer findMin(UniqueArray<T> array);
  void rotate3DVector(Particle particle, Real sin_Theta, Real cos_Theta, Real sin_Phi, Real cos_Phi);

public:

  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_TRACKINGMC(TrackingMCModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
