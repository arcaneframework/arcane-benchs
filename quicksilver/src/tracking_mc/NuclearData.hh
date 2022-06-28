// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef NUCLEAR_DATA_ARC_HH
#define NUCLEAR_DATA_ARC_HH

#include "arcane/ItemVector.h"
#include "structEnum.hh"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

using namespace Arcane;

class Polynomial
{
 public:
  Polynomial(Real aa, Real bb, Real cc, Real dd, Real ee)
  : _aa(aa)
  , _bb(bb)
  , _cc(cc)
  , _dd(dd)
  , _ee(ee)
  {}

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

  Real getCrossSection(Integer group);

  void sampleCollision(Real incidentEnergy, Real material_mass, RealUniqueArray& energyOut,
                       RealUniqueArray& angleOut, Integer& nOut, Int64* seed,
                       Integer max_production_size);

  RealUniqueArray _crossSection; //!< tabular data for microscopic cross section
  Enum _reactionType; //!< What type of reaction is this
  Real _nuBar; //!< If this is a fission, specify the nu bar
};

// This class holds an array of reactions for neutrons
class NuclearDataSpecies
{
 public:
  void addReaction(NuclearDataReaction::Enum type, Real nuBar,
                   RealUniqueArray energies, const Polynomial& polynomial,
                   Real reactionCrossSection);

  UniqueArray<NuclearDataReaction> _reactions;
};

// For this isotope, store the cross sections. In this case the species is just
// neutron.
class NuclearDataIsotope
{
 public:
  NuclearDataIsotope()
  : _species(1)
  {}

  UniqueArray<NuclearDataSpecies> _species;
};

// Top level class to handle all things related to nuclear data
class NuclearData
{
 public:
  NuclearData(Integer numGroups, Real energyLow, Real energyHigh);

  Integer addIsotope(Integer nReactions, const Polynomial& fissionFunction,
                     const Polynomial& scatterFunction,
                     const Polynomial& absorptionFunction, Real nuBar,
                     Real totalCrossSection, Real fissionWeight,
                     Real scatterWeight, Real absorptionWeight);

  Integer getEnergyGroup(Real energy);

  Integer getNumberReactions(Integer isotopeIndex);

  Real getTotalCrossSection(Integer isotopeIndex, Integer group);

  Real getReactionCrossSection(Integer reactIndex, Integer isotopeIndex,
                               Integer group);

  // Store the cross sections and reactions by isotope, which stores
  // it by species
  UniqueArray<NuclearDataIsotope> _isotopes;
  // This is the overall energy layout. If we had more than just
  // neutrons, this array would be a vector of vectors.
  RealUniqueArray _energies;
};

#endif