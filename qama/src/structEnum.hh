﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
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

// Si true, alors les résultats de sortie doivent être les mêmes que ceux de
// Quicksilver Original.
// (Mettre à true pour avoir le même comportement que QS Original).
#define QS_LEGACY_COMPATIBILITY true

// Permet d'activer le calcul des cross sections avant le tracking.
// (Mettre à false pour avoir le même comportement que QS Original).
#define QS_PRECOMPUTE_CROSS_SECTION false

// Permet d'activer la verification du statut des particules (pour prévenir le
// fait qu'une particule 'morte' soit toujours trackée, par exemple).
// Attention : Si l'option est désactivée, m_particle_status ne sera plus mis à
// 'ParticleState::oldParticle'. Il y aura donc une modification des valeurs et donc
// la comparaison bit-à-bit ne sera pas possible entre un write option true et un
// read option false, par exemple. Mis à part cette grandeur, il n'y aura pas
// de modification au niveau des résultats.
// (Mettre à false pour avoir le même comportement que QS Original).
#define QS_CHECK_PARTICLE_STATUS false

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

namespace PhysicalConstants
{

const Real _neutronRestMassEnergy = 9.395656981095e+2; /* MeV */
const Real _pi = 3.1415926535897932;
const Real _speedOfLight = 2.99792458e+10; // cm / s

// Constants used in math for computer science, roundoff, and other reasons
const Real _tinyDouble = 1.0e-13;
const Real _smallDouble = 1.0e-10;
const Real _hugeDouble = 1.0e+75;
//
} // namespace PhysicalConstants

#endif
