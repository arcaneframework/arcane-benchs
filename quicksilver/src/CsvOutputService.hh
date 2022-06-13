// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsvOutputService.hh                                         (C) 2000-2022 */
/*                                                                           */
/* Service permettant de construire et de sortir un .csv.                    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "ISimpleTableOutput.hh"
#include "arcane/BasicUnitTest.h"
#include "CsvOutput_axl.h"
#include <arcane/ServiceBuilder.h>
#include <arcane/ServiceFactory.h>

using namespace Arcane;

class CsvOutputService
: public ArcaneCsvOutputObject
{
public:
  CsvOutputService(const ServiceBuildInfo & sbi)
    : ArcaneCsvOutputObject(sbi)
    , m_size_rows(0)
    , m_size_columns(0)
    , m_name_rows(0)
    , m_name_columns(0)
    , m_path_computed(false)
    , m_name_tab_computed(false)
    , m_path_only_P0(true)
    , m_name_tab_only_P0(true)
    , m_last_row(-1)
    , m_last_column(-1)
    {
      if(sbi.creationType() == ST_CaseOption) {
        m_path = options()->getPath();
        m_with_option = true;
      }
      else {
        m_path = "./";
        m_with_option = false;
      }
      
    }
  
  virtual ~CsvOutputService() {};

public:
  void init() override;
  void init(String name_table) override;
  void init(String name_table, String separator) override;

  void clear() override;


  Integer addRow(String name_row) override;
  Integer addRow(String name_row, ConstArrayView<Real> elems) override;
  bool addRows(StringConstArrayView name_rows) override;

  Integer addColumn(String name_column) override;
  Integer addColumn(String name_column, ConstArrayView<Real> elems) override;
  bool addColumns(StringConstArrayView name_columns) override;
  

  bool addElemRow(Integer pos, Real elem) override;
  bool addElemRow(String name_row, Real elem, bool create_if_not_exist) override;
  bool addElemSameRow(Real elem) override;
  
  bool addElemsRow(Integer pos, ConstArrayView<Real> elems) override;
  bool addElemsRow(String name_row, ConstArrayView<Real> elems, bool create_if_not_exist) override;
  bool addElemsSameRow(ConstArrayView<Real> elems) override;


  bool addElemColumn(Integer pos, Real elem) override;
  bool addElemColumn(String name_column, Real elem, bool create_if_not_exist) override;
  bool addElemSameColumn(Real elem) override;

  bool addElemsColumn(Integer pos, ConstArrayView<Real> elems) override;
  bool addElemsColumn(String name_column, ConstArrayView<Real> elems, bool create_if_not_exist) override;
  bool addElemsSameColumn(ConstArrayView<Real> elems) override;


  bool editElemUp(Real elem) override;
  bool editElemDown(Real elem) override;
  bool editElemLeft(Real elem) override;
  bool editElemRight(Real elem) override;


  bool editElem(Integer pos_x, Integer pos_y, Real elem) override;
  bool editElem(String name_column, String name_row, Real elem) override;

  Real getElem(Integer pos_x, Integer pos_y) override;
  Real getElem(String name_column, String name_row) override;


  RealUniqueArray getRow(Integer pos) override;
  RealUniqueArray getColumn(Integer pos) override;
  
  RealUniqueArray getRow(String name_row) override;
  RealUniqueArray getColumn(String name_column) override;


  Integer getSizeRow(Integer pos) override;
  Integer getSizeColumn(Integer pos) override;
  
  Integer getSizeRow(String name_row) override;
  Integer getSizeColumn(String name_column) override;


  Integer getPosRow(String name_row) override;
  Integer getPosColumn(String name_column) override;

  Integer getNumRows() override;
  Integer getNumColumns() override;


  Integer addAverageColumn(String name_column) override;

  void print(bool only_P0) override;
  void print(Integer only_proc) override;
  bool writeFile(bool only_P0) override;
  bool writeFile(String path, bool only_P0) override;


private:
  String _computeAt(String name, bool& only_P0);
  String _computeFinal();
  void _print(std::ostream& stream);

private:
  String m_path;
  bool m_path_computed;
  bool m_path_only_P0;

  String m_name_tab;
  bool m_name_tab_computed;
  bool m_name_tab_only_P0;

  String m_separator;

  UniqueArray2<Real> m_values_csv;

  UniqueArray<String> m_name_rows;
  UniqueArray<String> m_name_columns;

  // Tailles des lignes/colonnes 
  // (et pas le nombre d'éléments, on compte les "trous" ici).
  UniqueArray<Integer> m_size_rows;
  UniqueArray<Integer> m_size_columns;

  // Dernier élement ajouté.
  Integer m_last_row;
  Integer m_last_column;

  bool m_with_option;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_CSVOUTPUT(CsvOutput, CsvOutputService);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
