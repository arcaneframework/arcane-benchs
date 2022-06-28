// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* IRandomNumberGenerator.hh                                   (C) 2000-2022 */
/*                                                                           */
/* Interface pour générateur de nombres aléatoires.                          */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef IRNG_HH
#define IRNG_HH

#include <arcane/ItemTypes.h>

using namespace Arcane;


class IRandomNumberGenerator
{
public:
  IRandomNumberGenerator() {};
  virtual ~IRandomNumberGenerator() {};
  
public:
  virtual void initSeed() = 0;
  virtual void initSeed(Int64 seed) = 0;

  virtual Int64 getSeed() = 0;

  virtual Int64 randomSeedGenerator() = 0;
  virtual Int64 randomSeedGenerator(Int64* parent_seed) = 0;

  virtual Real randomNumberGenerator() = 0;
  virtual Real randomNumberGenerator(Int64* seed) = 0;
};

#endif
