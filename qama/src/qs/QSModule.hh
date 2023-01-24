// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* QSModule.hh                                                 (C) 2000-2022 */
/*                                                                           */
/* Module principal QAMA                                                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/IMesh.h>
#include <arcane/IParallelMng.h>
#include <arcane/ITimeLoopMng.h>
#include <arcane/cartesianmesh/CellDirectionMng.h>
#include <arcane/cartesianmesh/ICartesianMesh.h>
#include <arcane/ICartesianMeshGenerationInfo.h>
#include <arcane/IMeshModifier.h>
#include <arcane/ServiceBuilder.h>
#include <arcane/ILoadBalanceMng.h>
#include <arcane/IMeshPartitionerBase.h>
#include <arcane/ISimpleTableOutput.h>
#include <arcane/ISimpleTableComparator.h>

#include "structEnum.hh"

enum eBoundaryCondition
{
  REFLECT,
  ESCAPE,
  OCTANT
};

#include "qs/QS_axl.h"


using namespace Arcane;

/*!
  \brief Module QS.
  Ce module permet d'initialiser les grandeurs au maillage et
  les tallies et permet d'afficher les résultats finaux.
 */
class QSModule : 
public ArcaneQSObject
{

 public:
  explicit QSModule(const ModuleBuildInfo& mbi)
  : ArcaneQSObject(mbi)
  {}

 public:
  void initModule() override;
  void startLoadBalancing() override;
  void cycleFinalize() override;
  void loopLoadBalancing() override;
  void afterLoadBalancing() override;
  void compareWithReference() override;
  void endModule() override;

  VersionInfo versionInfo() const override { return VersionInfo(1, 5, 0); }

 protected:
  void initMesh();
  void loadBalancing();
  ParticleEvent getBoundaryCondition(const Integer& pos);

protected:
  ISimpleTableOutput* m_csv;
  ISimpleTableComparator* m_csv_compare;
  Integer m_max_iteration;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_QS(QSModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
