// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "InitMC_axl.h"
#include <arcane/ITimeLoopMng.h>
#include <arcane/geometry/IGeometryMng.h>
#include <arcane/cartesianmesh/CellDirectionMng.h>
#include <arcane/cartesianmesh/ICartesianMesh.h>

using namespace Arcane;

/*!
 * \brief Module InitMC.
 */
class InitMCModule : public ArcaneInitMCObject {
public:
  explicit InitMCModule(const ModuleBuildInfo &mbi) : ArcaneInitMCObject(mbi) {}

public:

  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_INITMC(InitMCModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
