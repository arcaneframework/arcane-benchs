// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* structEnum.hh                                               (C) 2000-2022 */
/*                                                                           */
/* Struct/Enum en commun QAMA                                                */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef STRUCTENUM_HH
#define STRUCTENUM_HH

#include <arcane/ITimeLoopMng.h>

using namespace Arcane;

#define QS_LEGACY_COMPATIBILITY

enum eBoundaryCondition
{
  REFLECT,
  ESCAPE,
  OCTANT
};
enum eShape
{
  UNDEFINED,
  BRICK,
  SPHERE
};

enum CosDir
{
  MD_DirA = 0, // Alpha
  MD_DirB, // Beta
  MD_DirG // Gamma
};

enum ParticleEvent
{
  collision,
  census,
  faceEventUndefined,
  escape,
  reflection,
  cellChange,
  subDChange,
  undefined
};

enum ParticleState
{
  oldParticle, // Particle ayant déjà fait un tracking.
  newParticle, // Nouvelle particule.
  clonedParticle, // Particule provenant d'un split.
  exitedParticle, // Particule sortie du maillage.
  censusParticle // Particule ayant fini le tracking actuel.
};

#endif
