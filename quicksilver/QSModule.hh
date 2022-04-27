// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include <arcane/ITimeLoopMng.h>
#include <arcane/geometry/IGeometryMng.h>
#include <arcane/cartesianmesh/CellDirectionMng.h>
#include <arcane/cartesianmesh/ICartesianMesh.h>
#include <arcane/IParallelMng.h>
#include "arcane/IMesh.h"
#include "arcane/IMeshModifier.h"
#include "arcane/IItemFamily.h"

enum eBoundaryCondition{reflect, escape, octant}; // TODO : A deplacer (doit être defini avant QS_axl.h !).

#include "QS_axl.h"

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

/*!
 * \brief Module QS.
 */
class QSModule : 
public ArcaneQSObject {

public:
  explicit QSModule(const ModuleBuildInfo &mbi) : 
  ArcaneQSObject(mbi)

{
}

public:
  void initModule() override;
  void cycleFinalize() override;
  void endModule() override;

  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

public:

protected:
  ICartesianMesh* m_cartesian_mesh;


public:
  void cycleFinalizeTallies();
  void initMesh();
  void initTallies();
  UniqueArray<Face_Adjacency_Event> getBoundaryCondition();
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_QS(QSModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
