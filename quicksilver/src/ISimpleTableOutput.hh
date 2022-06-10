// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ISimpleTableOutput.hh                                       (C) 2000-2022 */
/*                                                                           */
/* Interface pour simples services de sortie de tableaux de valeurs.         */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef ISIMPLEOUTPUT_HH
#define ISIMPLEOUTPUT_HH

#include <arcane/ItemTypes.h>
#include <optional>

using namespace Arcane;


class ISimpleTableOutput
{
public:
  ISimpleTableOutput() {};
  virtual ~ISimpleTableOutput() {};
  
public:
  /**
   * @brief Méthode permettant d'initialiser le tableau.
   * 
   */
  virtual void init() = 0;
  
  /**
   * @brief Méthode permettant d'initialiser le tableau.
   * 
   * @param name_table Le nom du tableau (et du fichier de sortie).
   */
  virtual void init(String name_table) = 0;

  /**
   * @brief Méthode permettant d'initialiser le tableau.
   * 
   * @param name_table Le nom du tableau (et du fichier de sortie).
   * @param separator La chaine de caractère à inserer entre les valeurs lors de la sortie.
   */
  virtual void init(String name_table, String separator) = 0;

  /**
   * @brief Méthode permettant de nettoyer le tableau pour le réutiliser après.
   * 
   */
  virtual void clear() = 0;

  /**
   * @brief Méthode permettant d'ajouter une ligne.
   * 
   * @param name_row Le nom de la ligne.
   * @return Integer La position de la ligne dans le tableau.
   */
  virtual Integer addRow(String name_row) = 0;

  /**
   * @brief Méthode permettant d'ajouter une ligne.
   * 
   * @param name_row Le nom de la ligne.
   * @param elems Les éléments à insérer sur la ligne.
   * @return Integer La position de la ligne dans le tableau.
   */
  virtual Integer addRow(String name_row, ConstArrayView<Real> elems) = 0;

  /**
   * @brief Méthode permettant d'ajouter plusieurs lignes.
   * 
   * @param name_rows Les noms des lignes.
   * @return true Si toutes les lignes ont été créées.
   * @return false Si toutes les lignes n'ont pas été créées.
   */
  virtual bool addRows(StringConstArrayView name_rows) = 0;

  /**
   * @brief Méthode permettant d'ajouter une colonne.
   * 
   * @param name_column Le nom de la colonne.
   * @return Integer La position de la colonne dans le tableau.
   */
  virtual Integer addColumn(String name_column) = 0;

  /**
   * @brief Méthode permettant d'ajouter une colonne.
   * 
   * @param name_column Le nom de la colonne.
   * @param elems Les éléments à ajouter sur la colonne.
   * @return Integer La position de la colonne dans le tableau.
   */
  virtual Integer addColumn(String name_column, ConstArrayView<Real> elems) = 0;

  /**
   * @brief Méthode permettant d'ajouter plusieurs colonnes.
   * 
   * @param name_rows Les noms des colonnes.
   * @return true Si toutes les colonnes ont été créées.
   * @return false Si toutes les colonnes n'ont pas été créées.
   */
  virtual bool addColumns(StringConstArrayView name_columns) = 0;


  // Globalement, return true si tous les éléments ont été ajoutés.
  /**
   * @brief Méthode permettant d'ajouter un élément à une ligne.
   * 
   * @param pos La position de la ligne.
   * @param elem L'élément à ajouter.
   * @return true Si l'élément a pu être ajouté.
   * @return false Si l'élément n'a pas pu être ajouté.
   */
  virtual bool addElemRow(Integer pos, Real elem) = 0;
  /**
   * @brief Méthode permettant l'ajouter un élément sur une ligne.
   * 
   * @param name_row Le nom de la ligne.
   * @param elem L'élément à ajouter.
   * @param create_if_not_exist Pour savoir si l'on doit créer la ligne si elle n'existe pas encore.
   * @return true Si l'élément a pu être ajouté.
   * @return false Si l'élément n'a pas pu être ajouté.
   */
  virtual bool addElemRow(String name_row, Real elem, bool create_if_not_exist = true) = 0;
  /**
   * @brief Méthode permettant d'ajouter un élément sur la ligne dernièrement manipulée.
   * 
   * @param elem L'élément à ajouter.
   * @return true Si l'élément a été ajouté.
   * @return false Si l'élément n'a pas pu être ajouté.
   */
  virtual bool addElemSameRow(Real elem) = 0;

  /**
   * @brief Méthode permettant d'ajouter plusieurs éléments sur une ligne.
   * 
   * @param pos La position de la ligne.
   * @param elems Le tableau d'élement à ajouter.
   * @return true Si tous les éléments ont été ajoutés.
   * @return false Si [0;len(elems)[ éléments ont été ajoutés.
   */
  virtual bool addElemsRow(Integer pos, ConstArrayView<Real> elems) = 0;
  /**
   * @brief Méthode permettant d'ajouter plusieurs éléments sur une ligne.
   * 
   * @param name_row Le nom de la ligne.
   * @param elems Le tableau d'élement à ajouter.
   * @param create_if_not_exist Pour savoir si l'on doit créer la ligne si elle n'existe pas encore.
   * @return true Si tous les éléments ont été ajoutés.
   * @return false Si [0;len(elems)[ éléments ont été ajoutés.
   */
  virtual bool addElemsRow(String name_row, ConstArrayView<Real> elems, bool create_if_not_exist = true) = 0;
  /**
   * @brief Méthode permettant d'ajouter plusieurs éléments sur la ligne dernièrement manipulée.
   * 
   * @param elems Le tableau d'élement à ajouter.
   * @return true Si tous les éléments ont été ajoutés.
   * @return false Si [0;len(elems)[ éléments ont été ajoutés.
   */
  virtual bool addElemsSameRow(ConstArrayView<Real> elems) = 0;

  /**
   * @brief Méthode permettant d'ajouter un élément à une colonne.
   * 
   * @param pos La position de la colonne.
   * @param elem L'élément à ajouter.
   * @return true Si l'élément a pu être ajouté.
   * @return false Si l'élément n'a pas pu être ajouté.
   */
  virtual bool addElemColumn(Integer pos, Real elem) = 0;
  /**
   * @brief Méthode permettant l'ajouter un élément sur une colonne.
   * 
   * @param name_column Le nom de la colonne.
   * @param elem L'élément à ajouter.
   * @param create_if_not_exist Pour savoir si l'on doit créer la colonne si elle n'existe pas encore.
   * @return true Si l'élément a pu être ajouté.
   * @return false Si l'élément n'a pas pu être ajouté.
   */
  virtual bool addElemColumn(String name_column, Real elem, bool create_if_not_exist = true) = 0;
  /**
   * @brief Méthode permettant d'ajouter un élément sur la colonne dernièrement manipulée.
   * 
   * @param elem L'élément à ajouter.
   * @return true Si l'élément a été ajouté.
   * @return false Si l'élément n'a pas pu être ajouté.
   */
  virtual bool addElemSameColumn(Real elem) = 0;

  /**
   * @brief Méthode permettant d'ajouter plusieurs éléments sur une colonne.
   * 
   * @param pos La position de la colonne.
   * @param elems Le tableau d'élement à ajouter.
   * @return true Si tous les éléments ont été ajoutés.
   * @return false Si [0;len(elems)[ éléments ont été ajoutés.
   */
  virtual bool addElemsColumn(Integer pos, ConstArrayView<Real> elems) = 0;
  /**
   * @brief Méthode permettant d'ajouter plusieurs éléments sur une colonne.
   * 
   * @param name_column Le nom de la colonne.
   * @param elems Le tableau d'élement à ajouter.
   * @param create_if_not_exist Pour savoir si l'on doit créer la colonne si elle n'existe pas encore.
   * @return true Si tous les éléments ont été ajoutés.
   * @return false Si [0;len(elems)[ éléments ont été ajoutés.
   */
  virtual bool addElemsColumn(String name_column, ConstArrayView<Real> elems, bool create_if_not_exist = true) = 0;
  /**
   * @brief Méthode permettant d'ajouter plusieurs éléments sur la colonne dernièrement manipulée.
   * 
   * @param elems Le tableau d'élement à ajouter.
   * @return true Si tous les éléments ont été ajoutés.
   * @return false Si [0;len(elems)[ éléments ont été ajoutés.
   */
  virtual bool addElemsSameColumn(ConstArrayView<Real> elems) = 0;


  /**
   * @brief Méthode permettant de modifier un élement dans le tableau.
   * 
   * @param pos_x La position de la colonne à modifier.
   * @param pos_y La position de la ligne à modifier.
   * @param elem L'élement remplaçant.
   * @return true Si l'élement a bien été remplacé.
   * @return false Si l'élement n'a pas été remplacé.
   */
  virtual bool editElem(Integer pos_x, Integer pos_y, Real elem) = 0;
  /**
   * @brief Méthode permettant de modifier un élement dans le tableau.
   * 
   * @param name_column Le nom de la colonne où se trouve l'élement.
   * @param name_row Le nom de la ligne où se trouve l'élement.
   * @param elem L'élement remplaçant.
   * @return true Si l'élement a bien été remplacé.
   * @return false Si l'élement n'a pas pu être remplacé.
   */
  virtual bool editElem(String name_column, String name_row, Real elem) = 0;

  /**
   * @brief Méthode permettant d'obtenir une copie d'un élement.
   * 
   * @param pos_x La position de la colonne où se trouve l'élement.
   * @param pos_y La position de la ligne où se trouve l'élement.
   * @return Real L'élement trouvé (0 si non trouvé).
   */
  virtual Real getElem(Integer pos_x, Integer pos_y) = 0;
  /**
   * @brief Méthode permettant d'obtenir une copie d'un élement.
   * 
   * @param name_column Le nom de la colonne où se trouve l'élement.
   * @param name_row Le nom de la ligne où se trouve l'élement.
   * @return Real L'élement trouvé (0 si non trouvé).
   */
  virtual Real getElem(String name_column, String name_row) = 0;

  /**
   * @brief Méthode permettant d'obtenir une copie d'une ligne.
   * 
   * @param pos La position de la ligne.
   * @return RealUniqueArray La copie de la ligne (tableau vide si non trouvé).
   */
  virtual RealUniqueArray getRow(Integer pos) = 0;
  /**
   * @brief Méthode permettant d'obtenir une copie d'une ligne.
   * 
   * @param name_row Le nom de la ligne.
   * @return RealUniqueArray La copie de la ligne (tableau vide si non trouvé).
   */
  virtual RealUniqueArray getRow(String name_row) = 0;

  /**
   * @brief Méthode permettant d'obtenir une copie d'une colonne.
   * 
   * @param pos La position de la colonne.
   * @return RealUniqueArray La copie de la colonne (tableau vide si non trouvé).
   */
  virtual RealUniqueArray getColumn(Integer pos) = 0;
  /**
   * @brief Méthode permettant d'obtenir une copie d'une colonne.
   * 
   * @param name_column Le nom de la colonne.
   * @return RealUniqueArray La copie de la colonne (tableau vide si non trouvé).
   */
  virtual RealUniqueArray getColumn(String name_column) = 0;

  /**
   * @brief Méthode permettant d'obtenir la taille d'une ligne.
   * Incluant les hypotétiques trous dans la ligne.
   * 
   * @param pos La position de la ligne.
   * @return Integer La taille de la ligne (0 si non trouvé).
   */
  virtual Integer getSizeRow(Integer pos) = 0;
  /**
   * @brief Méthode permettant d'obtenir la taille d'une ligne.
   * Incluant les hypotétiques trous dans la ligne.
   * 
   * @param pos Le nom de la ligne.
   * @return Integer La taille de la ligne (0 si non trouvé).
   */
  virtual Integer getSizeRow(String name_row) = 0;
  
  /**
   * @brief Méthode permettant d'obtenir la taille d'une colonne.
   * Incluant les hypotétiques trous dans la colonne.
   * 
   * @param pos La position de la colonne.
   * @return Integer La taille de la colonne (0 si non trouvé).
   */
  virtual Integer getSizeColumn(Integer pos) = 0;
  /**
   * @brief Méthode permettant d'obtenir la taille d'une colonne.
   * Incluant les hypotétiques trous dans la colonne.
   * 
   * @param pos Le nom de la colonne.
   * @return Integer La taille de la colonne (0 si non trouvé).
   */
  virtual Integer getSizeColumn(String name_column) = 0;

  /**
   * @brief Méthode permettant de récupérer la position d'une ligne.
   * 
   * @param name_row Le nom de la ligne.
   * @return Integer La position de la ligne (-1 si non trouvé).
   */
  virtual Integer getPosRow(String name_row) = 0;
  /**
   * @brief Méthode permettant de récupérer la position d'une colonne.
   * 
   * @param name_row Le nom de la colonne.
   * @return Integer La position de la colonne (-1 si non trouvé).
   */
  virtual Integer getPosColumn(String name_column) = 0;
  
  /**
   * @brief Méthode permettant de récupérer le nombre de lignes dans le tableau.
   * C'est, en quelque sorte, le max d'élements que peut contenir une colonne.
   * 
   * @return Integer Le nombre de lignes du tableau.
   */
  virtual Integer getNumRows() = 0;
  /**
   * @brief Méthode permettant de récupérer le nombre de colonnes dans le tableau.
   * C'est, en quelque sorte, le max d'élements que peut contenir une ligne.
   * 
   * @return Integer Le nombre de colonnes du tableau.
   */
  virtual Integer getNumColumns() = 0;
  
  /**
   * @brief Méthode permettant de créer une colonne contenant la moyenne des éléments de chaque ligne
   * 
   * @param name_column Le nom de la nouvelle colonne.
   * @return Integer La position de la colonne.
   */
  virtual Integer addAverageColumn(String name_column) = 0;

  /**
   * @brief Méthode permettant d'afficher le tableau.
   * 
   * @param only_P0 S'il y a uniquement le P0 qui doit l'afficher.
   */
  virtual void print(bool only_P0 = true) = 0;
  /**
   * @brief Méthode permettant d'afficher le tableau.
   * 
   * @param only_proc L'id du processus devant afficher le tableau.
   */
  virtual void print(Integer only_proc) = 0;
  /**
   * @brief Méthode permettant d'écrire le tableau dans un fichier.
   * Si only_P0 == true, return indéfini.
   * 
   * @param only_P0 S'il y a uniquement le P0 qui doit l'écrire.
   * @return true Si le fichier a été correctement écrit.
   * @return false Si le fichier n'a pas été correctement écrit.
   */
  virtual bool writeFile(bool only_P0 = false) = 0;
  /**
   * @brief Méthode permettant d'écrire le tableau dans un fichier.
   * Si only_P0 == true, return indéfini.
   * 
   * @param path Le répertoire où sera écrit le fichier.
   * @param only_P0 S'il y a uniquement le P0 qui doit l'écrire.
   * @return true Si le fichier a été correctement écrit.
   * @return false Si le fichier n'a pas été correctement écrit.
   */
  virtual bool writeFile(String path, bool only_P0 = false) = 0;
};

#endif