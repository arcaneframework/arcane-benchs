#include "NuclearData.hh"
#include <cmath>
#include "MC_RNG_State.hh"
#include "qs_assert.hh"

using std::log10;
using std::pow;

// Set the cross section values and reaction type
// Cross sections are scaled to produce the supplied reactionCrossSection at 1MeV.
NuclearDataReaction::NuclearDataReaction(
   Enum reactionType, Real nuBar, const RealUniqueArray energies,
   const Polynomial& polynomial, Real reactionCrossSection)
: _crossSection(energies.size()-1),
  _reactionType(reactionType),
  _nuBar(nuBar)
{
   Integer nGroups = _crossSection.size();

   for (Integer ii=0; ii<nGroups; ++ii)
   {
      Real energy = (energies[ii] + energies[ii+1]) / 2.0;
      _crossSection[ii] = pow( 10, polynomial(log10( energy)));
   }

   // Find the normalization value for the polynomial.  This is the
   // value of the energy group that contains 1 MeV
   Real normalization = 1.0;
   for (Integer ii=0; ii<nGroups; ++ii)
      if (energies[ii+1] > 1. ) //1 MeV
      {
         normalization = _crossSection[ii];
         break;
      }
   qs_assert(normalization > 0.);

   // scale to specified reaction cross section
   Real scale = reactionCrossSection/normalization;
   for (Integer ii=0; ii<nGroups; ++ii)
      _crossSection[ii] *= scale;
}



void NuclearDataReaction::sampleCollision(
   Real incidentEnergy, Real material_mass, Real* energyOut,
   Real* angleOut, Integer &nOut, Int64* seed, Integer max_production_size)
{
   Real randomNumber;
   switch(_reactionType)
   {
     case Scatter:
      nOut = 1;
      randomNumber = rngSample(seed);
      energyOut[0] = incidentEnergy * (1.0 - (randomNumber*(1.0/material_mass)));
      randomNumber = rngSample(seed) * 2.0 - 1.0;
      angleOut[0] = randomNumber;
      break;
     case Absorption:
      break;
     case Fission:
      {
         Integer numParticleOut = (Integer)(_nuBar + rngSample(seed));
         qs_assert( numParticleOut <= max_production_size );
         nOut = numParticleOut;
         for (Integer outIndex = 0; outIndex < numParticleOut; outIndex++)
         {
            randomNumber = rngSample(seed) / 2.0 + 0.5;
            energyOut[outIndex] = (20 * randomNumber*randomNumber);
            randomNumber = rngSample(seed) * 2.0 - 1.0;
            angleOut[outIndex] = randomNumber;
         }
      }
      break;
     case Undefined:
      printf("_reactionType invalid\n");
      qs_assert(false);
   }
}



// Then call this for each reaction to set cross section values
void NuclearDataSpecies::addReaction(
   NuclearDataReaction::Enum type, Real nuBar,
   RealUniqueArray energies, const Polynomial& polynomial, Real reactionCrossSection)
{
   _reactions.add(NuclearDataReaction(type, nuBar, energies, polynomial, reactionCrossSection));
}



// Set up the energies boundaries of the neutron
NuclearData::NuclearData(Integer numGroups, Real energyLow, Real energyHigh) : _energies( numGroups+1)
{
   qs_assert (energyLow < energyHigh);
   _energies[0] = energyLow;
   _energies[numGroups] = energyHigh;
   Real logLow = log(energyLow);
   Real logHigh = log(energyHigh);
   Real delta = (logHigh - logLow) / (numGroups + 1.0);
   for (Integer energyIndex = 1; energyIndex < numGroups; energyIndex++)
   {
      Real logValue = logLow + delta *energyIndex;
      _energies[energyIndex] = exp(logValue);
   }
}

Integer NuclearData::addIsotope(
   Integer nReactions,
   const Polynomial& fissionFunction,
   const Polynomial& scatterFunction,
   const Polynomial& absorptionFunction,
   Real nuBar,
   Real totalCrossSection,
   Real fissionWeight, Real scatterWeight, Real absorptionWeight)
{
   _isotopes.add(NuclearDataIsotope());

   Real totalWeight = fissionWeight + scatterWeight + absorptionWeight;

   Integer nFission    = nReactions / 3;
   Integer nScatter    = nReactions / 3;
   Integer nAbsorption = nReactions / 3;
   switch (nReactions % 3)
   {
     case 0:
      break;
     case 1:
      ++nScatter;
      break;
     case 2:
      ++nScatter;
      ++nFission;
      break;
   }
   
   Real fissionCrossSection    = (totalCrossSection * fissionWeight)    / (nFission    * totalWeight);
   Real scatterCrossSection    = (totalCrossSection * scatterWeight)    / (nScatter    * totalWeight);
   Real absorptionCrossSection = (totalCrossSection * absorptionWeight) / (nAbsorption * totalWeight);

   _isotopes.back()._species[0]._reactions.reserve( nReactions);

   for (Integer ii=0; ii<nReactions; ++ii)
   {
      NuclearDataReaction::Enum type;
      Polynomial polynomial(0.0, 0.0, 0.0, 0.0, 0.0);
      Real reactionCrossSection = 0.;
      // reaction index % 3 is one of the 3 reaction types
      switch (ii % 3)
      {
        case 0:
         type = NuclearDataReaction::Scatter;
         polynomial = scatterFunction;
         reactionCrossSection = scatterCrossSection;
         break;
        case 1:
         type = NuclearDataReaction::Fission;
         polynomial = fissionFunction;
         reactionCrossSection = fissionCrossSection;
         break;
        case 2:
         type = NuclearDataReaction::Absorption;
         polynomial = absorptionFunction;
         reactionCrossSection = absorptionCrossSection;
         break;
      }
      _isotopes.back()._species[0].addReaction(type, nuBar, _energies, polynomial, reactionCrossSection);
   }
   

   return _isotopes.size() - 1;
}


// Return the cross section for this energy group
Real NuclearDataReaction::getCrossSection(Integer group)
{
   qs_assert(group < _crossSection.size());
   return _crossSection[group];
}



Integer NuclearData::getNumberReactions(Integer isotopeIndex)
{
   qs_assert(isotopeIndex < _isotopes.size());
   return _isotopes[isotopeIndex]._species[0]._reactions.size();
}


// For this energy, return the group index

Integer NuclearData::getEnergyGroup(Real energy)
{
   Integer numEnergies = _energies.size();
   if (energy <= _energies[0]) return 0;
   if (energy > _energies[numEnergies-1]) return numEnergies-1;

   Integer high = numEnergies-1;
   Integer low = 0;

   while( high != low+1 )
   {
       Integer mid = (high+low)/2;
       if( energy < _energies[mid] ) 
           high = mid;
       else
           low  = mid;
   }

   return low;
}


// General routines to help access data lower down
// Return the total cross section for this energy group

Real NuclearData::getTotalCrossSection(Integer isotopeIndex, Integer group)
{
   qs_assert(isotopeIndex < _isotopes.size());
   Integer numReacts = _isotopes[isotopeIndex]._species[0]._reactions.size();
   Real totalCrossSection = 0.0;
   for (Integer reactIndex = 0; reactIndex < numReacts; reactIndex++)
   {
      totalCrossSection += _isotopes[isotopeIndex]._species[0]._reactions[reactIndex].getCrossSection(group);
   }
   return totalCrossSection;
}


// Return the total cross section for this energy group

Real NuclearData::getReactionCrossSection(
   Integer reactIndex, Integer isotopeIndex, Integer group)
{
   qs_assert(isotopeIndex < _isotopes.size());
   qs_assert(reactIndex < _isotopes[isotopeIndex]._species[0]._reactions.size());
   return _isotopes[isotopeIndex]._species[0]._reactions[reactIndex].getCrossSection(group);
}


