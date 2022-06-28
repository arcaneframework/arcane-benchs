// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* IRandomNumberGenerator.hh                                   (C) 2000-2022 */
/*                                                                           */
/* Interface pour générateur de nombres aléatoires.                          */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef IRNG_HH
#define IRNG_HH

#include <arcane/ItemTypes.h>

using namespace Arcane;


class IRandomNumberGenerator
{
public:
  IRandomNumberGenerator() {};
  virtual ~IRandomNumberGenerator() {};
  
public:
  /**
   * @brief Méthode permettant d'initialiser le service.
   * 
   * Avec la seed en option (ou la seed par défaut si l'on est
   * en mode singleton).
   * 
   */
  virtual void initSeed() = 0;

  /**
   * @brief Méthode permettant d'initialiser le service.
   * 
   * @param seed La seed d'origine.
   */
  virtual void initSeed(Int64 seed) = 0;

  /**
   * @brief Méthode permettant de récupérer la seed actuelle.
   * 
   * @return Int64 La seed.
   */
  virtual Int64 getSeed() = 0;

  /**
   * @brief Méthode permettant de générer une autre seed à partir de
   * la seed en mémoire.
   * 
   * @return Int64 La nouvelle seed.
   */
  virtual Int64 randomSeedGenerator() = 0;

  /**
   * @brief Méthode permettant de générer une autre seed à partir de
   * la seed transmise en paramètre.
   * 
   * Cette méthode n'utilise pas la seed en mémoire.
   * 
   * @param parent_seed La seed d'origine.
   * @return Int64 La nouvelle seed.
   */
  virtual Int64 randomSeedGenerator(Int64* parent_seed) = 0;

  /**
   * @brief Méthode permettant de générer un nombre aléatoire avec
   * la seed en mémoire.
   * 
   * @return Real Le nombre généré.
   */
  virtual Real randomNumberGenerator() = 0;

   /**
   * @brief Méthode permettant de générer un nombre aléatoire avec
   * la seed transmise en paramètre.
   *    
   * Cette méthode n'utilise pas la seed en mémoire.
   * 
   * @param seed La seed.
   * @return Real Le nombre généré.
   */
  virtual Real randomNumberGenerator(Int64* seed) = 0;
};

#endif
