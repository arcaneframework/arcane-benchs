// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsvOutputService.cc                                         (C) 2000-2022 */
/*                                                                           */
/* Service permettant de construire et de sortir un .csv.                    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "CsvOutputService.hh"
using namespace Arcane;

void CsvOutputService::
init(String name_csv, String separator_csv)
{
  separator = separator_csv;
  rows[0] = name_csv + separator;
  size_rows.add(1);
  size_columns.add(1);
}

Integer CsvOutputService::
addRow(String name_row, bool fill_start)
{
  name_rows.add(name_row);
  String new_line = name_row + separator;

  if(fill_start) {
    for(Integer i = 1; i < size_rows[0]; i++)
      new_line = new_line + separator;
  }

  size_rows.add(size_columns[0]);

  rows.add(new_line);
  return size_columns[0]++;
}

Integer CsvOutputService::
addRow(String name_row, ConstArrayView<Real>& elems)
{
  name_rows.add(name_row);
  String new_line = name_row + separator;

  size_rows.add(0);

  rows.add(new_line);
  addElemsRow(size_columns[0], elems);
  return size_columns[0]++;
}

bool CsvOutputService::
addElemsRow(Integer pos, ConstArrayView<Real>& elems)
{
  if(elems.size() + size_rows[pos] > size_rows[0]) warning() << "Attention, tous les élements ne seront pas mis.";

  String real_string;

  for(Integer i = 0; i < size_rows[0] && i < elems.size(); i++) {
    real_string = real_string + String::fromNumber(elems[i]) + separator;
    size_rows[pos]++;
  }
  rows[pos] = rows[pos] + real_string;
  return true;
}

bool CsvOutputService::
addElemRow(Integer pos, Real elem)
{
  if(pos >= size_columns[0]) {
    error() << "Mauvaise pos";
    return false;
  }

  rows[pos] = rows[pos] + String::fromNumber(elem) + separator;
  size_rows[pos]++;
  return true;
}

bool CsvOutputService::
addElemRow(String name_row, Real elem, bool create_if_not_exist)
{
  std::optional<Integer> pos = name_rows.span().findFirst(name_row);

  if(pos)                       return addElemRow(pos.value()+1, elem);
  else if(create_if_not_exist)  return addElemRow(addRow(name_row, false), elem);
  else                          return false;
}

Integer CsvOutputService::
addColumn(String name_column, bool fill_start)
{
  name_columns.add(name_column);
  String new_column = name_column + separator;

  size_columns.add(1);

  rows[0] = rows[0] + new_column;
  return size_rows[0]++;
}

Integer CsvOutputService::
addColumn(String name_column, ConstArrayView<Real>& elems)
{
  name_columns.add(name_column);
  String new_column = name_column + separator;

  size_columns.add(1);

  rows[0] = rows[0] + new_column;
  addElemsColumn(size_rows[0], elems);
  return size_rows[0]++;
}

bool CsvOutputService::
addElemsColumn(Integer pos, ConstArrayView<Real>& elems)
{
  if(elems.size() + size_columns[pos] > size_columns[0]) warning() << "Attention, tous les élements ne seront pas mis.";

  for(Integer i = size_columns[pos], j = 0; i < size_columns[0] && j < elems.size(); i++, j++) {
    rows[i] = rows[i] + String::fromNumber(elems[j]) + separator;
    size_columns[pos]++;
  }
  return true;
}

void CsvOutputService::
print()
{
  for(Integer i = 0; i < rows.size(); i++) {
    std::cout << rows[i] << std::endl;
  }
}

void CsvOutputService::
writeFile(String name_file)
{
  std::ofstream ofile( (options()->getPath() + name_file).localstr() );
  for(Integer i = 0; i < rows.size(); i++) {
    ofile << rows[i] << std::endl;
  }
  ofile.close();
}
