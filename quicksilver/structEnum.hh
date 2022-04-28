#ifndef STRUCTENUM_HH
#define STRUCTENUM_HH

using namespace Arcane;


enum eBoundaryCondition{reflect, escape, octant};
enum eShape{UNDEFINED, BRICK, SPHERE};

enum cosDir
{
  MD_DirA = 0, // Alpha
  MD_DirB,     // Beta
  MD_DirG      // Gamma
};

enum faceEvent
{
  Collision,
  Facet_Crossing_Transit_Exit,
  Census,
  Facet_Crossing_Tracking_Error,
  Facet_Crossing_Escape,
  Facet_Crossing_Reflection,
  Facet_Crossing_Communication
};

enum faceAdjacencyEvent
{
  Adjacency_Undefined = 0,
  Boundary_Escape,
  Boundary_Reflection,
  Transit_On_Processor,
  Transit_Off_Processor
};

enum segmentOutcomeType
{
  Initialize                    = -1,
  Collision1                     = 0,
  Facet_Crossing                = 1,
  Census1                        = 2,
  Max_Number                    = 3
};

enum particleEvent
{
  undefined,
  collision,
  census,
  faceEscape,
  faceReflection,
  cellChange,
  sdChange,
  errorEvent
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

static bool ordre_qs[]  = {true, true, false, false, false, true};

static Integer QS2ArcaneFace[] = {4, 1, 5, 2, 3, 0};

static Integer QS2ArcaneNode[] = {0, 1, 2, 3,
                       0, 3, 2, 1,
                       1, 0, 3, 2,
                       0, 1, 2, 3,
                       0, 1, 2, 3,
                       0, 3, 2, 1};


#endif
