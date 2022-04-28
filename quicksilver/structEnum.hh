#ifndef STRUCTENUM_HH
#define STRUCTENUM_HH

using namespace Arcane;


enum eBoundaryCondition{REFLECT, ESCAPE, OCTANT};
enum eShape{UNDEFINED, BRICK, SPHERE};

enum cosDir
{
  MD_DirA = 0, // Alpha
  MD_DirB,     // Beta
  MD_DirG      // Gamma
};

enum particleEvent
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

struct Nearest_Facet
{
   Integer facet;
   Real distance_to_facet;
   Real dot_product;
   
   Nearest_Facet()
   : facet(0),
     distance_to_facet(1e80),
     dot_product(0.0)
   {}
};

struct Distance_To_Facet
{
    Real distance;
    Integer facet;
    Integer subfacet;
    
    Distance_To_Facet()
    : distance(0.0),
      facet(0),
      subfacet(0) 
    {}
};

static bool scan_order[]  = {true, true, false, false, false, true};

static Integer QS2ArcaneFace[] = {4, 1, 5, 2, 3, 0};

static Integer QS2ArcaneNode[] = {0, 1, 2, 3,
                       0, 3, 2, 1,
                       1, 0, 3, 2,
                       0, 1, 2, 3,
                       0, 1, 2, 3,
                       0, 3, 2, 1};


#endif
