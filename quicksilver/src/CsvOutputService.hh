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
    , rows(1)
    , size_rows(0)
    , size_columns(0)
    , name_rows(0)
    , name_columns(0)
    {}
  
  virtual ~CsvOutputService() {};

public:
  virtual void init(String name_csv, String separator);

  virtual Integer addRow(String name_row, bool fill_start);
  virtual Integer addRow(String name_row, ConstArrayView<Real>& elems);
  virtual Integer addColumn(String name_column, bool fill_start);
  virtual Integer addColumn(String name_column, ConstArrayView<Real>& elems);
  
  virtual bool addElemRow(Integer pos, Real elem);
  virtual bool addElemRow(String name_row, Real elem, bool create_if_not_exist);

  virtual bool addElemColumn(Integer pos, Real elem){return false;}
  virtual bool addElemColumn(String name_column, Real elem, bool create_if_not_exist){return false;}

  virtual void print();
  virtual void writeFile(String name_file);

private:
  bool addElemsRow(Integer pos, ConstArrayView<Real>& elems);
  bool addElemsColumn(Integer pos, ConstArrayView<Real>& elems);

private:
  UniqueArray<String> rows;

  UniqueArray<Integer> size_rows;
  UniqueArray<Integer> size_columns;

  UniqueArray<String> name_rows;     // TODO : Mettre HashTableMapT
  UniqueArray<String> name_columns;  // TODO : Mettre HashTableMapT

  String separator;
};

#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_CSVOUTPUT(CsvOutput, CsvOutputService);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/