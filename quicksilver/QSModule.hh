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
#include "QS_axl.h"



#include "CoralBenchmark.hh"
#include "CycleTracking.hh"
#include "EnergySpectrum.hh"
#include "MC_Fast_Timer.hh"
#include "MC_Particle_Buffer.hh"
#include "MC_Processor_Info.hh"
#include "MC_SourceNow.hh"
#include "MC_Time_Info.hh"
#include "MonteCarlo.hh"
#include "NVTX_Range.hh"
#include "Parameters.hh"
#include "ParticleVault.hh"
#include "ParticleVaultContainer.hh"
#include "PopulationControl.hh"
#include "SendQueue.hh"
#include "Tallies.hh"
#include "cudaFunctions.hh"
#include "cudaUtils.hh"
#include "initMC.hh"
#include "macros.hh"
#include "qs_assert.hh"
#include "utils.hh"
#include "utilsMpi.hh"
#include "MC_Segment_Outcome.hh"
#include "MC_Nearest_Facet.hh"
#include "MC_Distance_To_Facet.hh"



#include "git_hash.hh"
#include "git_vers.hh"

using namespace Arcane;
using namespace std;

/*!
 * \brief Module QS.
 */
class QSModule : 
public ArcaneQSObject {
public:
  explicit QSModule(const ModuleBuildInfo &mbi) : 
  ArcaneQSObject(mbi)
, m_mesh(mbi.mesh())
, m_particle_family(nullptr)
, m_particle_family_with_ghost(nullptr)
, m_first_uid(0) 
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

  IMesh* m_mesh;
  IItemFamily* m_particle_family;
  IItemFamily* m_particle_family_with_ghost;

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

  Int64 m_first_uid;
  SharedArray< SharedArray<Integer> > m_extra_ghost_particles_to_send;

  ParticleVectorView m_processingView;

  MonteCarlo *monteCarlo = NULL;
  MonteCarlo *monteCarloArc = NULL;
  Parameters params;

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
  ICartesianMesh* cartesian_mesh;

public:
  void CycleFinalizeTallies();
  void getParametersAxl();
  MonteCarlo* initMCArc(const Parameters& params);
  void initNuclearData(MonteCarlo* monteCarlo, const Parameters& params);
  void initMesh(MonteCarlo* monteCarlo, const Parameters& params);
  qs_vector<MC_Subfacet_Adjacency_Event::Enum> getBoundaryCondition(const Parameters& params);
  void initTallies(MonteCarlo* monteCarlo, const Parameters& params);
  void initializeCentersRandomly(int nCenters,
                                const GlobalFccGrid& grid,
                                vector<MC_Vector>& centers);
  void initializeCentersGrid(double lx, double ly, double lz,
                            int xDom, int yDom, int zDom,
                            vector<MC_Vector>& centers);
  void consistencyCheck(int myRank, const qs_vector<MC_Domain>& domain);
  bool isInside(const GeometryParameters& geom, const MC_Vector& rr);
  string findMaterial(const Parameters& params, const MC_Vector& rr);
  void checkCrossSections(MonteCarlo* monteCarlo, const Parameters& params);
  void clearCrossSectionCache();
  void MC_SourceNowArc(MonteCarlo *monteCarlo);
  Real Get_Speed_From_Energy(Particle p);
  void MCT_Generate_Coordinate_3D_GArc(Particle p,
                                  // uint64_t *random_number_seed,
                                  // Cell &cell,
                                  // MC_Vector &coordinate,  
                                  MonteCarlo* monteCarlo );
  double MCT_Cell_Volume_3D_G_vector_tetDetArc(const MC_Vector &v0_,
                                            const MC_Vector &v1_,
                                            const MC_Vector &v2_,
                                            const MC_Vector &v3);

  void PopulationControlArc(MonteCarlo* monteCarlo, bool loadBalance);
  void PopulationControlGutsArc(const double splitRRFactor, uint64_t currentNumParticles);
  void RouletteLowWeightParticlesArc(MonteCarlo* monteCarlo);

  void tracking(MonteCarlo* monteCarlo);
  void trackingArc(MonteCarlo* monteCarlo);
  void CollisionEventSuite();
  void CycleTrackingGutsArc( MonteCarlo *monteCarlo, Particle particle );
  void CycleTrackingFunctionArc( MonteCarlo *monteCarlo, Particle particle);
  MC_Segment_Outcome_type::Enum MC_Segment_OutcomeArc(MonteCarlo* monteCarlo, Particle particle, unsigned int &flux_tally_index);
  double weightedMacroscopicCrossSectionArc(MonteCarlo* monteCarlo, Cell cell, int energyGroup);
  double macroscopicCrossSectionArc(MonteCarlo* monteCarlo, int reactionIndex, Cell cell, int isoIndex, int energyGroup);
  MC_Nearest_Facet MCT_Nearest_FacetArc(Particle particle,
                                      double distance_threshold,
                                      double current_best_distance,
                                      bool new_segment,
                                      MonteCarlo* monteCarlo );
  MC_Nearest_Facet MCT_Nearest_Facet_3D_GArc( Particle particle);
  double MCT_Nearest_Facet_3D_G_Distance_To_Segment(double plane_tolerance,
                                                     double facet_normal_dot_direction_cosine,
                                                     double A, double B, double C, double D,
                                                     const MC_Vector &facet_coords0,
                                                     const MC_Vector &facet_coords1,
                                                     const MC_Vector &facet_coords2,
                                                     Particle particle,
                                                     bool allow_enter);

  MC_Nearest_Facet MCT_Nearest_Facet_Find_NearestArc(Particle particle,
                                  int &iteration, // input/output
                                  double &move_factor, // input/output
                                  int num_facets_per_cell,
                                  MC_Distance_To_Facet *distance_to_facet,
                                  int &retry /* output */ );
  MC_Nearest_Facet MCT_Nearest_Facet_Find_Nearest(int num_facets_per_cell,
                                                   MC_Distance_To_Facet *distance_to_facet);

  void MCT_Nearest_Facet_3D_G_Move_ParticleArc(Particle particle, // input/output: move this coordinate
                                          double move_factor);
  int CollisionEventArc(MonteCarlo* monteCarlo, Particle particle);
  void updateTrajectory( double energy, double angle, Particle particle );
  MC_Tally_Event::Enum MC_Facet_Crossing_EventArc(Particle particle, MonteCarlo* monteCarlo);
  void MCT_Reflect_ParticleArc(MonteCarlo *monteCarlo, Particle particle);
  unsigned int MC_Find_Min(const double *array, int num_elements);
  void Sample_Isotropic(Particle p);
  void copyParticle(Particle pSrc, Particle pNew);
  void copyParticles(Int32UniqueArray idsSrc, Int32UniqueArray idsNew);
  void Rotate3DVector(Particle particle, double sin_Theta, double cos_Theta, double sin_Phi, double cos_Phi);
  void initParticle(Particle p, int64_t rns);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_QS(QSModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
