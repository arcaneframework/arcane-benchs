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

#ifndef CSVOUTPUTSERVICE_HH
#define CSVOUTPUTSERVICE_HH

#include "ISimpleOutput.hh"
#include "CsvOutput_axl.h"

using namespace Arcane;

class CsvOutputService
: public ArcaneCsvOutputObject
{
public:
  CsvOutputService(const ServiceBuildInfo & sbi)
    : ArcaneCsvOutputObject(sbi) 
    , m_rows(1)
    , m_size_rows(0)
    , m_size_columns(0)
    , m_name_rows(0)
    , m_name_columns(0)
    , m_name_computed(false)
    , m_only_P0(true)
    {
      if(sbi.creationType() == ST_CaseOption) {
        m_path_name = options()->getFile();
      }
      else {
        m_path_name = "./output.csv";
      }
    }
  
  virtual ~CsvOutputService() {};

public:
  void init(String name_csv, String separator) override;

  Integer addRow(String name_row, bool fill_start) override;
  Integer addRow(String name_row, ConstArrayView<Real>& elems) override;
  Integer addColumn(String name_column, bool fill_start) override;
  Integer addColumn(String name_column, ConstArrayView<Real>& elems) override;
  
  bool addElemRow(Integer pos, Real elem) override;
  bool addElemRow(String name_row, Real elem, bool create_if_not_exist) override;

  bool addElemColumn(Integer pos, Real elem) override {return false;}
  bool addElemColumn(String name_column, Real elem, bool create_if_not_exist) override {return false;}

  void print(bool only_P0) override;
  bool writeFile() override;
  bool writeFile(String path_file) override;

private:
  bool addElemsRow(Integer pos, ConstArrayView<Real>& elems);
  bool addElemsColumn(Integer pos, ConstArrayView<Real>& elems);
  void computePathName();

private:
  UniqueArray<String> m_rows;

  UniqueArray<Integer> m_size_rows;
  UniqueArray<Integer> m_size_columns;

  UniqueArray<String> m_name_rows;
  UniqueArray<String> m_name_columns;

  String m_separator;

  String m_path_name;
  bool m_name_computed;
  bool m_only_P0;
};

#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_CSVOUTPUT(CsvOutput, CsvOutputService);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/