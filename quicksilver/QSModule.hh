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

#include "QS_axl.h"


#include "NuclearDataArc.hh"
#include "MC_Vector.hh"



using namespace Arcane;

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
 * \brief Module QS.
 */
class QSModule : 
public ArcaneQSObject {

public:
  explicit QSModule(const ModuleBuildInfo &mbi) : 
  ArcaneQSObject(mbi)
, m_particle_family(nullptr)
, m_cycle(0)

{
}

public:
  void startInit() override;

  void cycleInit() override;
  void cycleTracking() override;
  void cycleFinalize() override;

  void gameOver() override;

  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

public:

  IItemFamily* m_particle_family;

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

  ParticleVectorView m_processingView;

  Int64 m_nx, m_ny, m_nz;

  std::atomic<Int64> m_absorb_a{0};
  std::atomic<Int64> m_census_a{0};
  std::atomic<Int64> m_escape_a{0};
  std::atomic<Int64> m_collision_a{0};
  std::atomic<Int64> m_end_a{0};
  std::atomic<Int64> m_fission_a{0};
  std::atomic<Int64> m_produce_a{0};
  std::atomic<Int64> m_scatter_a{0};
  std::atomic<Int64> m_start_a{0};
  std::atomic<Int64> m_source_a{0};
  std::atomic<Int64> m_rr_a{0};
  std::atomic<Int64> m_split_a{0};
  std::atomic<Int64> m_numSegments_a{0};

protected:
  ICartesianMesh* m_cartesian_mesh;
  Arcane::Materials::IMeshMaterialMng* material_mng;
  NuclearDataArc* m_nuclearData;

  Real m_source_particle_weight;
  Integer m_cycle; // TODO : Voir pour utiliser var arcane interne.

public:
  void CycleFinalizeTallies();
  void getParametersAxl();
  void initMCArc();
  void initNuclearDataArc();
  void initMesh();
  UniqueArray<Face_Adjacency_Event> getBoundaryCondition();
  void initTallies();
  void clearCrossSectionCache();
  void MC_SourceNowArc();
  Real Get_Speed_From_Energy(Particle p);
  void MCT_Generate_Coordinate_3D_GArc(Particle p);
  Real MCT_Cell_Volume_3D_G_vector_tetDetArc(const MC_Vector &v0_,
                                            const MC_Vector &v1_,
                                            const MC_Vector &v2_,
                                            const MC_Vector &v3);

  void PopulationControlArc();
  void PopulationControlGutsArc(const Real splitRRFactor, Int64 currentNumParticles);
  void RouletteLowWeightParticlesArc();

  void trackingArc();
  void CollisionEventSuite();
  void CycleTrackingGutsArc( Particle particle );
  void CycleTrackingFunctionArc( Particle particle);
  Segment_Outcome_type MC_Segment_OutcomeArc(Particle particle);
  Real weightedMacroscopicCrossSectionArc(Cell cell, int energyGroup);
  Real macroscopicCrossSectionArc(int reactionIndex, Cell cell, int isoIndex, int energyGroup);
  Nearest_Facet MCT_Nearest_FacetArc(Particle particle,
                                      Real distance_threshold,
                                      Real current_best_distance,
                                      bool new_segment);
  Nearest_Facet MCT_Nearest_Facet_3D_GArc( Particle particle);
  Real MCT_Nearest_Facet_3D_G_Distance_To_Segment(Real plane_tolerance,
                                                     Real facet_normal_dot_direction_cosine,
                                                     Real A, Real B, Real C, Real D,
                                                     const MC_Vector &facet_coords0,
                                                     const MC_Vector &facet_coords1,
                                                     const MC_Vector &facet_coords2,
                                                     Particle particle,
                                                     bool allow_enter);

  Nearest_Facet MCT_Nearest_Facet_Find_NearestArc(Particle particle,
                                  int &iteration, // input/output
                                  Real &move_factor, // input/output
                                  int num_facets_per_cell,
                                  Distance_To_Facet *distance_to_facet,
                                  int &retry /* output */ );
  Nearest_Facet MCT_Nearest_Facet_Find_Nearest(int num_facets_per_cell,
                                                   Distance_To_Facet *distance_to_facet);

  void MCT_Nearest_Facet_3D_G_Move_ParticleArc(Particle particle, // input/output: move this coordinate
                                          Real move_factor);
  int CollisionEventArc(Particle particle);
  void updateTrajectory( Real energy, Real angle, Particle particle );
  Tally_Event MC_Facet_Crossing_EventArc(Particle particle);
  void MCT_Reflect_ParticleArc(Particle particle);
  Integer MC_Find_Min(RealUniqueArray array);
  void Sample_Isotropic(Particle p);
  void copyParticle(Particle pSrc, Particle pNew);
  void copyParticles(Int32UniqueArray idsSrc, Int32UniqueArray idsNew);
  void Rotate3DVector(Particle particle, Real sin_Theta, Real cos_Theta, Real sin_Phi, Real cos_Phi);
  void initParticle(Particle p, Int64 rns);
  bool isInsideArc(Integer posOptions, Cell cell);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_QS(QSModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
