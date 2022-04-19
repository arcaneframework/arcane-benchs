// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include <arcane/IParticleExchanger.h>
#include <arcane/ITimeLoopMng.h>
#include <arcane/geometry/IGeometryMng.h>
#include <arcane/cartesianmesh/CellDirectionMng.h>
#include <arcane/cartesianmesh/ICartesianMesh.h>
#include <arcane/IParallelMng.h>
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
class QSModule : public ArcaneQSObject {
public:
  explicit QSModule(const ModuleBuildInfo &mbi) : ArcaneQSObject(mbi) {}

public:
  void startInit() override;

  void cycleInit() override;
  void cycleTracking() override;
  void cycleFinalize() override;

  void gameOver() override;

  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

public:
  MonteCarlo *monteCarlo = NULL;
  MonteCarlo *monteCarloArc = NULL;
  Parameters params;

protected:
  ICartesianMesh* cartesian_mesh;

public:
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
  double Get_Speed_From_Energy(double energy);
  void MCT_Generate_Coordinate_3D_GArc(uint64_t *random_number_seed,
                                  Cell &cell,
                                  MC_Vector &coordinate,  
                                  MonteCarlo* monteCarlo );
  double MCT_Cell_Volume_3D_G_vector_tetDetArc(const MC_Vector &v0_,
                                            const MC_Vector &v1_,
                                            const MC_Vector &v2_,
                                            const MC_Vector &v3);

  void PopulationControlArc(MonteCarlo* monteCarlo, bool loadBalance);
  void PopulationControlGutsArc(const double splitRRFactor, uint64_t currentNumParticles, ParticleVaultContainer* my_particle_vault);
  void RouletteLowWeightParticlesArc(MonteCarlo* monteCarlo);

  void tracking(MonteCarlo* monteCarlo);
  void trackingArc(MonteCarlo* monteCarlo);
  void CycleTrackingGutsArc( MonteCarlo *monteCarlo, int particle_index, ParticleVault *processingVault, ParticleVault *processedVault );
  void CycleTrackingFunctionArc( MonteCarlo *monteCarlo, MC_Particle &mc_particle, int particle_index, ParticleVault* processingVault, ParticleVault* processedVault);
  MC_Segment_Outcome_type::Enum MC_Segment_OutcomeArc(MonteCarlo* monteCarlo, MC_Particle &mc_particle, unsigned int &flux_tally_index);
  double weightedMacroscopicCrossSectionArc(MonteCarlo* monteCarlo, Cell &cell, int energyGroup);
  double macroscopicCrossSectionArc(MonteCarlo* monteCarlo, int reactionIndex, Cell& cell, int isoIndex, int energyGroup);
  MC_Nearest_Facet MCT_Nearest_FacetArc(MC_Particle *mc_particle,
                                      MC_Vector &coordinate,
                                      const DirectionCosine *direction_cosine,
                                      double distance_threshold,
                                      double current_best_distance,
                                      bool new_segment,
                                      MonteCarlo* monteCarlo );
  MC_Nearest_Facet MCT_Nearest_Facet_3D_GArc( MC_Particle *mc_particle,
                                      MC_Vector &coordinate,
                                      const DirectionCosine *direction_cosine);
  double MCT_Nearest_Facet_3D_G_Distance_To_Segment(double plane_tolerance,
                                                     double facet_normal_dot_direction_cosine,
                                                     double A, double B, double C, double D,
                                                     const MC_Vector &facet_coords0,
                                                     const MC_Vector &facet_coords1,
                                                     const MC_Vector &facet_coords2,
                                                     const MC_Vector &coordinate,
                                                     const DirectionCosine *direction_cosine,
                                                     bool allow_enter);

  MC_Nearest_Facet MCT_Nearest_Facet_Find_NearestArc(MC_Particle *mc_particle,
                                  MC_Vector &coordinate,
                                  int &iteration, // input/output
                                  double &move_factor, // input/output
                                  int num_facets_per_cell,
                                  MC_Distance_To_Facet *distance_to_facet,
                                  int &retry /* output */ );
  MC_Nearest_Facet MCT_Nearest_Facet_Find_Nearest(int num_facets_per_cell,
                                                   MC_Distance_To_Facet *distance_to_facet);

  void MCT_Nearest_Facet_3D_G_Move_ParticleArc(Cell& cell,
                                          MC_Vector &coordinate, // input/output: move this coordinate
                                          double move_factor);
  bool CollisionEventArc(MonteCarlo* monteCarlo, MC_Particle &mc_particle);
  void updateTrajectory( double energy, double angle, MC_Particle& particle );
  MC_Tally_Event::Enum MC_Facet_Crossing_EventArc(MC_Particle &mc_particle, MonteCarlo* monteCarlo, int particle_index, ParticleVault* processingVault);
  void MCT_Reflect_ParticleArc(MonteCarlo *monteCarlo, MC_Particle &particle);
  unsigned int MC_Find_Min(const double *array, int num_elements);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_QS(QSModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
