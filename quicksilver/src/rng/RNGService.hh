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

#include "IRandomNumberGenerator.hh"
#include "rng/RNG_axl.h"

using namespace Arcane;

class RNGService
: public ArcaneRNGObject
{
public:
  RNGService(const ServiceBuildInfo & sbi)
    : ArcaneRNGObject(sbi)
    , m_seed(0)
    {
      m_with_option = (sbi.creationType() == ST_CaseOption);
    }
  
  virtual ~RNGService() {};

public:
  void initSeed() override;
  void initSeed(Int64 seed) override;

  Int64 getSeed() override;

  Int64 randomSeedGenerator() override;
  Int64 randomSeedGenerator(Int64* parent_seed) override;

  Real randomNumberGenerator() override;
  Real randomNumberGenerator(Int64* seed) override;

protected:
  void _breakupUInt64(uint64_t uint64_in, uint32_t& front_bits, uint32_t& back_bits);
  uint64_t _reconstructUInt64(uint32_t front_bits, uint32_t back_bits);
  void _pseudoDES(uint32_t& lword, uint32_t& irword);
  uint64_t _hashState(uint64_t initial_number);

protected:
  Int64 m_seed;
  bool m_with_option;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_RNG(RNG, RNGService);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
