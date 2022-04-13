// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "QS_axl.h"
#include <arcane/ITimeLoopMng.h>
#include <arcane/geometry/IGeometryMng.h>
#include <arcane/cartesianmesh/CellDirectionMng.h>
#include <arcane/cartesianmesh/ICartesianMesh.h>


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
  void checkCrossSections(MonteCarlo* monteCarlo, const Parameters& params);

};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_QS(QSModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
