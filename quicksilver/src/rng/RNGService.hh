// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* RNGService.hh                                               (C) 2000-2022 */
/*                                                                           */
/* Implémentation d'un générateur de nombres aléatoires.                     */
/* Basé sur le générateur de Quicksilver (LLNL).                             */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/IRandomNumberGenerator.h>
#include "rng/RNG_axl.h"

using namespace Arcane;

class RNGService
: public ArcaneRNGObject
{
public:
  RNGService(const ServiceBuildInfo & sbi)
    : ArcaneRNGObject(sbi)
    , m_seed(0)
    , m_size_of_seed(sizeof(Int64))
    {
    }
  
  virtual ~RNGService() {};

public:
  bool initSeed() override;
  bool initSeed(RandomNumberGeneratorSeed seed) override;

  RandomNumberGeneratorSeed seed() override;
  RandomNumberGeneratorSeed emptySeed() override;

  Integer neededSizeOfSeed() override;

  bool isLeapSeedSupported() override { return false; };
  RandomNumberGeneratorSeed generateRandomSeed(Integer leap) override;
  RandomNumberGeneratorSeed generateRandomSeed(RandomNumberGeneratorSeed* parent_seed, Integer leap) override;

  bool isLeapNumberSupported() override { return false; };
  Real generateRandomNumber(Integer leap) override;
  Real generateRandomNumber(RandomNumberGeneratorSeed* seed, Integer leap) override;

protected:
  Real _rngSample(Int64* seed);
  void _breakupUInt64(uint64_t uint64_in, uint32_t& front_bits, uint32_t& back_bits);
  uint64_t _reconstructUInt64(uint32_t front_bits, uint32_t back_bits);
  void _pseudoDES(uint32_t& lword, uint32_t& irword);
  uint64_t _hashState(uint64_t initial_number);

protected:
  Int64 m_seed;
  Integer m_size_of_seed;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_RNG(RNG, RNGService);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
