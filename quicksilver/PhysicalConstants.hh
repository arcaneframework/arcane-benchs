#ifndef PHYSICAL_CONSTANTS_HH
#define PHYSICAL_CONSTANTS_HH

#include "DeclareMacro.hh"
HOST_DEVICE_CLASS
namespace PhysicalConstants
{

const Real _neutronRestMassEnergy = 9.395656981095e+2; /* MeV */
const Real _pi = 3.1415926535897932;
const Real _speedOfLight  = 2.99792458e+10;                // cm / s

// Constants used in math for computer science, roundoff, and other reasons
 const Real _tinyDouble           = 1.0e-13;
 const Real _smallDouble          = 1.0e-10;
 const Real _hugeDouble           = 1.0e+75;
//
}
HOST_DEVICE_END


#endif
