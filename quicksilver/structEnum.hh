#ifndef STRUCTENUM_HH
#define STRUCTENUM_HH

#include <arcane/utils/Real3.h>

using namespace Arcane;

#define MAX_PRODUCTION_SIZE 4

enum eBoundaryCondition{REFLECT, ESCAPE, OCTANT};
enum eShape{UNDEFINED, BRICK, SPHERE};

enum CosDir
{
  MD_DirA = 0, // Alpha
  MD_DirB,     // Beta
  MD_DirG      // Gamma
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
  oldParticle,    // Particle ayant déjà fait un tracking.
  newParticle,    // Nouvelle particule.
  clonedParticle, // Particule provenant d'un split.
  exitedParticle, // Particule sortie du maillage.
  censusParticle  // Particule ayant fini le tracking actuel.
};

struct DistanceToFacet
{
    Real distance;
    Integer facet;
    
    DistanceToFacet()
    : distance(0.0),
      facet(0)
    {}
};

static bool scan_order[]  = {true, true, false, false, false, true};

static Integer QS_to_arcaneFace[] = {4, 1, 5, 2, 3, 0};

static Integer QS_to_arcaneNode[] = { 0, 1, 2, 3,
                                      0, 3, 2, 1,
                                      1, 0, 3, 2,
                                      0, 1, 2, 3,
                                      0, 1, 2, 3,
                                      0, 3, 2, 1};

static Real3 avToReal3(RealArrayView av)
{
  return Real3(av[MD_DirX], av[MD_DirY],av[MD_DirZ]);
}

#endif
