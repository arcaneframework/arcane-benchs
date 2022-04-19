#ifndef MC_RNG_STATE_INCLUDE
#define MC_RNG_STATE_INCLUDE

#include "portability.hh"
#include "DeclareMacro.hh"
#include <math.h>

//----------------------------------------------------------------------------------------------------------------------
//  A random number generator that implements a 64 bit linear congruential generator (lcg).
//
//  This implementation is based on the rng class from Nick Gentile.
//----------------------------------------------------------------------------------------------------------------------

// Generate a new random number seed
HOST_DEVICE
int64_t rngSpawn_Random_Number_Seed(int64_t *parent_seed);
uint64_t rngSpawn_Random_Number_Seed(uint64_t *parent_seed);
HOST_DEVICE_END

//----------------------------------------------------------------------------------------------------------------------
//  Sample returns the pseudo-random number produced by a call to a random
//  number generator.
//----------------------------------------------------------------------------------------------------------------------

HOST_DEVICE
inline double rngSample(int64_t *seed)
{
   // Reset the state from the previous value.
   *seed = 2862933555777941757L*(*seed) + 3037000493L;
   *seed &= ~(1UL << 63);

   // Map the int state in (0,2**64) to double (0,1)
   // by multiplying by
   // 1/(2**64 - 1) = 1/18446744073709551615.
   return 5.4210108624275222e-20*(*seed);
}
HOST_DEVICE_END

HOST_DEVICE
inline double rngSample(uint64_t *seed)
{
   // Reset the state from the previous value.
   *seed = 2862933555777941757L*(*seed) + 3037000493L;
   *seed &= ~(1UL << 63);

   // Map the int state in (0,2**64) to double (0,1)
   // by multiplying by
   // 1/(2**64 - 1) = 1/18446744073709551615.
   return 5.4210108624275222e-20*(*seed);
}
HOST_DEVICE_END

#endif
