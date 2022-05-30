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
#include <arcane/ServiceBuilder.h>
#include "ISimpleOutput.hh"

#include "structEnum.hh"

#include "QS_axl.h"

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
  void cycleFinalize() override;
  void endModule() override;

  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

 protected:
  void initMesh();
  ParticleEvent getBoundaryCondition(const Integer& pos);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_QS(QSModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
