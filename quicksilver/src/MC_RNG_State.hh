// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef MC_RNG_STATE_INCLUDE
#define MC_RNG_STATE_INCLUDE

#ifdef CSTDINT_MISSING
#include <stdint.h>
#else
#include <cstdint>
#endif

#include <arcane/ITimeLoopMng.h>
#include <math.h>

using namespace Arcane;

//----------------------------------------------------------------------------------------------------------------------
//  A random number generator that implements a 64 bit linear congruential
//  generator (lcg).
//
//  This implementation is based on the rng class from Nick Gentile.
//----------------------------------------------------------------------------------------------------------------------

// Generate a new random number seed

Int64 rngSpawn_Random_Number_Seed(Int64* parent_seed);

//----------------------------------------------------------------------------------------------------------------------
//  Sample returns the pseudo-random number produced by a call to a random
//  number generator.
//----------------------------------------------------------------------------------------------------------------------

inline Real
rngSample(Int64* seed)
{
  // Reset the state from the previous value.
  *seed = 2862933555777941757ULL * (uint64_t)(*seed) + 3037000493ULL;
  //*seed &= ~(1UL << 63);
  // Map the int state in (0,2**64) to double (0,1)
  // by multiplying by
  // 1/(2**64 - 1) = 1/18446744073709551615.
  volatile Real fin = 5.4210108624275222e-20 * (uint64_t)(*seed);
  ARCANE_ASSERT(fin >= 0, "rngSample negative");
  return fin;
}

#endif
