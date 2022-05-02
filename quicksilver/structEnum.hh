#ifndef STRUCTENUM_HH
#define STRUCTENUM_HH

using namespace Arcane;

#define MAX_PRODUCTION_SIZE 4
#define LOG false

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

struct NearestFacet
{
   Integer facet;
   Real distance_to_facet;
   Real dot_product;
   
   NearestFacet()
   : facet(0),
     distance_to_facet(1e80),
     dot_product(0.0)
   {}
};

struct DistanceToFacet
{
    Real distance;
    Integer facet;
    Integer subfacet;
    
    DistanceToFacet()
    : distance(0.0),
      facet(0),
      subfacet(0) 
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


#endif
