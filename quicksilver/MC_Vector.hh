#ifndef MC_VECTOR_INCLUDE
#define MC_VECTOR_INCLUDE

#include "arcane/ItemVector.h"
#include <arcane/utils/Real3.h>
#include <cmath>

using namespace Arcane;

class MC_Vector
: public Real3
{
public:
  MC_Vector(Real a, Real b, Real c)
  {
    x = a;
    y = b;
    z = c;
  }

  MC_Vector()
  {
    x = 0;
    y = 0;
    z = 0;
  }

  MC_Vector(RealArrayView av)
  {
    x = av[MD_DirX];
    y = av[MD_DirY];
    z = av[MD_DirZ];
  }

  MC_Vector(const Real3& av)
  {
    x = av.x;
    y = av.y;
    z = av.z;
  }

  // Distance from this vector to another point.
  inline double
  Distance(const MC_Vector& vv) const
  {
    Real3 fin = *this - vv;
    return fin.normL2();
  }

  inline double
  Dot(const MC_Vector& tmp) const
  {
    Real3 fin = *this * tmp;
    return fin[0] + fin[1] + fin[2];
  }

  inline MC_Vector
  Cross(const MC_Vector& v) const
  {
    return MC_Vector(
      y * v.z - z * v.y,
      z * v.x - x * v.z,
      x * v.y - y * v.x);
  }
};

#endif
