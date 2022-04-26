#ifndef NUCLEAR_DATA_ARC_HH
#define NUCLEAR_DATA_ARC_HH

#include <cstdio>
#include <string>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "qs_assert.hh"
#include "DeclareMacro.hh"
#include "arcane/ItemVector.h"

using namespace Arcane;


class Polynomial
{
 public:
   Polynomial(Real aa, Real bb, Real cc, Real dd, Real ee)
   :
   _aa(aa), _bb(bb), _cc(cc), _dd(dd), _ee(ee){}

   Real operator()(Real xx) const
   {
      return _ee + xx * (_dd + xx * (_cc + xx * (_bb + xx * (_aa))));
   }

 private:
   Real _aa, _bb, _cc, _dd, _ee;
};

// Lowest level class at the reaction level
class NuclearDataReaction
{
 public:
   // The types of reactions
   enum Enum
   {
      Undefined = 0,
      Scatter,
      Absorption,
      Fission
   };
   
   NuclearDataReaction(){};

   NuclearDataReaction(Enum reactionType, Real nuBar, RealUniqueArray energies,
                       const Polynomial& polynomial, Real reationCrossSection);
   

   HOST_DEVICE_CUDA
   Real getCrossSection(Integer group);
   HOST_DEVICE_CUDA
   void sampleCollision(Real incidentEnergy, Real material_mass, Real* energyOut,
                        Real* angleOut, Integer &nOut, Int64* seed, Integer max_production_size);
   
   RealUniqueArray _crossSection; //!< tabular data for microscopic cross section
   Enum _reactionType;                //!< What type of reaction is this
   Real _nuBar;                     //!< If this is a fission, specify the nu bar

};

// This class holds an array of reactions for neutrons
class NuclearDataSpecies
{
 public:
   
   void addReaction(NuclearDataReaction::Enum type, Real nuBar, RealUniqueArray energies,
                    const Polynomial& polynomial, Real reactionCrossSection);
   
   UniqueArray<NuclearDataReaction> _reactions;
};

// For this isotope, store the cross sections. In this case the species is just neutron.
class NuclearDataIsotope
{
 public:
   
   NuclearDataIsotope()
   : _species(1){}
   
   UniqueArray<NuclearDataSpecies> _species;

};

// Top level class to handle all things related to nuclear data
class NuclearData
{
 public:
   
   NuclearData(Integer numGroups, Real energyLow, Real energyHigh);

   Integer addIsotope(Integer nReactions,
                  const Polynomial& fissionFunction,
                  const Polynomial& scatterFunction,
                  const Polynomial& absorptionFunction,
                  Real nuBar,
                  Real totalCrossSection,
                  Real fissionWeight, Real scatterWeight, Real absorptionWeight);

   HOST_DEVICE_CUDA
   Integer getEnergyGroup(Real energy);
   HOST_DEVICE_CUDA
   Integer getNumberReactions(Integer isotopeIndex);
   HOST_DEVICE_CUDA
   Real getTotalCrossSection(Integer isotopeIndex, Integer group);
   HOST_DEVICE_CUDA
   Real getReactionCrossSection(Integer reactIndex, Integer isotopeIndex, Integer group);

   // Store the cross sections and reactions by isotope, which stores
   // it by species
   UniqueArray<NuclearDataIsotope> _isotopes;
   // This is the overall energy layout. If we had more than just
   // neutrons, this array would be a vector of vectors.
   RealUniqueArray _energies;

};

#endif

// The input for the nuclear data comes from the material section
// The input looks may like
//
// material NAME
// nIsotope=XXX
// nReactions=XXX
// fissionCrossSection="XXX"
// scatterCrossSection="XXX"
// absorptionCrossSection="XXX"
// nuBar=XXX
// totalCrossSection=XXX
// fissionWeight=XXX
// scatterWeight=XXX
// absorptionWeight=XXX
//
// Material NAME2
// ...
//
// table NAME
// a=XXX
// b=XXX
// c=XXX
// d=XXX
// e=XXX
//
// table NAME2
//
// Each isotope inside a material will have identical cross sections.
// However, it will be treated as unique in the nuclear data.
// Cross sectionsare strings that refer to tables
