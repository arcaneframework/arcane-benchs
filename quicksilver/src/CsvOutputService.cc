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
#include <arcane/IMesh.h>
#include <arcane/IParallelMng.h>
#include <optional>

using namespace Arcane;

void CsvOutputService::
init(String name_csv, String separator_csv)
{
  m_separator = separator_csv;
  m_rows[0] = name_csv + m_separator;
  m_size_rows.add(1);
  m_size_columns.add(1);
}

Integer CsvOutputService::
addRow(String name_row, bool fill_start)
{
  m_name_rows.add(name_row);
  String new_line = name_row + m_separator;

  if(fill_start) {
    for(Integer i = 1; i < m_size_rows[0]; i++)
      new_line = new_line + m_separator;
    m_size_rows.add(m_size_columns[0]);
  }
  else {
    m_size_rows.add(0);
  }


  m_rows.add(new_line);
  return m_size_columns[0]++;
}

Integer CsvOutputService::
addRow(String name_row, ConstArrayView<Real>& elems)
{
  m_name_rows.add(name_row);
  String new_line = name_row + m_separator;

  m_size_rows.add(0);

  m_rows.add(new_line);
  addElemsRow(m_size_columns[0], elems);
  return m_size_columns[0]++;
}

bool CsvOutputService::
addElemsRow(Integer pos, ConstArrayView<Real>& elems)
{
  if(elems.size() + m_size_rows[pos] > m_size_rows[0]) warning() << "Attention, tous les élements ne seront pas mis.";

  String real_string;

  for(Integer i = 0; i < m_size_rows[0] && i < elems.size(); i++) {
    real_string = real_string + String::fromNumber(elems[i]) + m_separator;
    m_size_rows[pos]++;
  }
  m_rows[pos] = m_rows[pos] + real_string;
  return true;
}

bool CsvOutputService::
addElemRow(Integer pos, Real elem)
{
  if(pos >= m_size_columns[0]) {
    error() << "Mauvaise pos";
    return false;
  }

  m_rows[pos] = m_rows[pos] + String::fromNumber(elem) + m_separator;
  m_size_rows[pos]++;
  return true;
}

bool CsvOutputService::
addElemRow(String name_row, Real elem, bool create_if_not_exist)
{
  std::optional<Integer> pos = m_name_rows.span().findFirst(name_row);

  if(pos)                       return addElemRow(pos.value()+1, elem);
  else if(create_if_not_exist)  return addElemRow(addRow(name_row, false), elem);
  else                          return false;
}

Integer CsvOutputService::
addColumn(String name_column, bool fill_start)
{
  m_name_columns.add(name_column);
  String new_column = name_column + m_separator;

  m_size_columns.add(1);

  m_rows[0] = m_rows[0] + new_column;
  return m_size_rows[0]++;
}

Integer CsvOutputService::
addColumn(String name_column, ConstArrayView<Real>& elems)
{
  m_name_columns.add(name_column);
  String new_column = name_column + m_separator;

  m_size_columns.add(1);

  m_rows[0] = m_rows[0] + new_column;
  addElemsColumn(m_size_rows[0], elems);
  return m_size_rows[0]++;
}

bool CsvOutputService::
addElemsColumn(Integer pos, ConstArrayView<Real>& elems)
{
  if(elems.size() + m_size_columns[pos] > m_size_columns[0]) warning() << "Attention, tous les élements ne seront pas mis.";

  for(Integer i = m_size_columns[pos], j = 0; i < m_size_columns[0] && j < elems.size(); i++, j++) {
    m_rows[i] = m_rows[i] + String::fromNumber(elems[j]) + m_separator;
    m_size_columns[pos]++;
  }
  return true;
}

void CsvOutputService::
computePathName()
{
  // On ne lit le nom qu'une seule fois.
  if(m_name_computed) return;
  m_name_computed = true;

  // On découpe la string là où se trouve les @.
  StringUniqueArray string_splited;
  m_path_name.split(string_splited, '@');

  // On traite les mots entre les "@".
  if(string_splited.size() > 1) {
    // On recherche "proc_id" dans le tableau (donc @proc_id@ dans le nom).
    std::optional<Integer> proc_id = string_splited.span().findFirst("proc_id");
    // On remplace "@proc_id@" par l'id du proc.
    if(proc_id)
    {
      string_splited[proc_id.value()] = String::fromNumber(mesh()->parallelMng()->commRank());
      m_only_P0 = false;
    }
    // Il n'y a que P0 qui write.
    else{
      m_only_P0 = true;
    }
  }

  // Autres traitements.
  if(!m_path_name.endsWith(".csv")) {
    string_splited.add(".csv");
  }

  // On recombine la chaine.
  String combined = "";
  for(String str : string_splited){
    combined = combined + str;
  }
  m_path_name = combined;
}

void CsvOutputService::
print(bool only_P0)
{
  computePathName();
  if(only_P0 && mesh()->parallelMng()->commRank() != 0) return;
  pinfo() << "P" << mesh()->parallelMng()->commRank() << " - Ecriture du .csv dans la sortie standard :";
  for(Integer i = 0; i < m_rows.size(); i++) {
    std::cout << m_rows[i] << std::endl;
  }
  pinfo() << "P" << mesh()->parallelMng()->commRank() << " - Fin écriture .csv";
}

bool CsvOutputService::
writeFile()
{
  computePathName();
  if(m_only_P0 && mesh()->parallelMng()->commRank() != 0) return true;

  std::ofstream ofile(m_path_name.localstr());
  if(ofile.fail()) return false;
  for(Integer i = 0; i < m_rows.size(); i++) {
    ofile << m_rows[i] << std::endl;
  }
  ofile.close();
  return true;
}

bool CsvOutputService::
writeFile(String path_file)
{
  m_path_name = path_file;
  m_name_computed = false;
  return writeFile();
}