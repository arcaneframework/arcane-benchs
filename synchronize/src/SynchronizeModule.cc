// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* SynchronizeModule.cc                                         (C) 2000-2022 */
/*                                                                           */
/* Bench Synchronize.                                                         */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/ISubDomain.h"
#include "arcane/IMesh.h"
#include "arcane/MathUtils.h"
#include "arcane/ITimeLoopMng.h"
#include "arcane/VariableTypes.h"
#include "arcane/ItemEnumerator.h"
#include "arcane/IParallelMng.h"
#include "arcane/ModuleFactory.h"
#include "arcane/ItemPrinter.h"
#include "arcane/ITimeStats.h"
#include "arcane/accelerator/core/IAcceleratorMng.h"

#include "arcane/mesh/ItemFamily.h"

#include "arcane/UnstructuredMeshConnectivity.h"

#include "arcane/accelerator/Reduce.h"
#include "arcane/accelerator/Runner.h"
#include "arcane/accelerator/VariableViews.h"
#include "arcane/accelerator/RunCommandEnumerate.h"

#include "Synchronize_axl.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Synchronize
{

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Module Synchronize.
 */
class SynchronizeModule
: public ArcaneSynchronizeObject
{
 public:

  //! Constructeur
  explicit SynchronizeModule(const ModuleBuildInfo& mb);

 public:

  VersionInfo versionInfo() const override { return VersionInfo(2, 0, 1); }

 public:

  void doOneIteration() override;

 private:
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

SynchronizeModule::
SynchronizeModule(const ModuleBuildInfo& sbi)
: ArcaneSynchronizeObject(sbi)
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SynchronizeModule::
doOneIteration()
{
  info() << "Launching infinite loop of synchronize...";
  VariableCellReal3 cell_test1(VariableBuildInfo(meshHandle(),"TMP1"));
  cell_test1.fill(Real3(1.0,1.0,1.0));

  VariableCellReal cell_test2(VariableBuildInfo(meshHandle(),"TMP2"));
  cell_test2.fill(1.0);

  VariableCellReal3 cell_test3(VariableBuildInfo(meshHandle(),"TMP3"));
  cell_test3.fill(Real3(1.0,1.0,1.0));

  VariableCellReal cell_test4(VariableBuildInfo(meshHandle(),"TMP4"));
  cell_test4.fill(1.0);

  VariableCellReal3 cell_test5(VariableBuildInfo(meshHandle(),"TMP5"));
  cell_test5.fill(Real3(1.0,1.0,1.0));

  Int64 nb_iteration = 0;
  while (1){
    ++nb_iteration;
    if ((nb_iteration%50000)==0)
      info() << "Iteration " << nb_iteration << " ...";
    cell_test1.synchronize();
    cell_test2.synchronize();
    cell_test3.synchronize();
    cell_test4.synchronize();
    cell_test5.synchronize();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_SYNCHRONIZE(SynchronizeModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // End namespace Synchronize

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
