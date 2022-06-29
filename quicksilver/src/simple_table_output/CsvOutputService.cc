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
#include <arcane/Directory.h>

#include <optional>


using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsvOutputService::
init()
{
  if(m_with_option && options()->getTableName() != "") {
    init(options()->getTableName());
  }
  else {
    init("Table_@proc_id@");
  }

}

void CsvOutputService::
init(String name_table)
{
  m_name_tab = _computeAt(name_table, m_name_tab_only_P0);
  m_name_tab_computed = true;

  m_separator = ";";

  if(m_with_option && options()->getTableDir() != "") {
    m_path = _computeAt(options()->getTableDir(), m_path_only_P0);
  }
  else {
    m_path = _computeAt("./csv/", m_path_only_P0);
  }
  m_path_computed = true;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsvOutputService::
clear()
{
  m_values_csv.clear();

  m_name_rows.clear();
  m_name_columns.clear();

  m_size_rows.clear();
  m_size_columns.clear();

  m_last_row = -1;
  m_last_column = -1;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Integer CsvOutputService::
addRow(String name_row)
{
  Integer pos = m_values_csv.dim1Size();
  m_values_csv.resize(pos+1);

  m_name_rows.add(name_row);
  m_size_rows.add(0);

  m_last_row = pos;

  return pos;
}

Integer CsvOutputService::
addRow(String name_row, ConstArrayView<Real> elems)
{
  Integer pos = m_values_csv.dim1Size();
  m_values_csv.resize(pos+1);
    
  m_name_rows.add(name_row);
  m_size_rows.add(0);

  addElemsRow(pos, elems);

  return pos;
}

bool CsvOutputService::
addRows(StringConstArrayView name_rows)
{
  Integer size = name_rows.size();
  if(size == 0) return true;

  Integer pos = m_values_csv.dim1Size();
  m_values_csv.resize(pos+size);

  m_name_rows.addRange(name_rows);
  m_size_rows.addRange(IntegerUniqueArray(size, 0));

  m_last_row = pos;

  return true;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Integer CsvOutputService::
addColumn(String name_column)
{
  Integer pos = m_values_csv.dim2Size();
  m_values_csv.resize(m_values_csv.dim1Size(), pos+1);

  m_name_columns.add(name_column);
  m_size_columns.add(0);

  m_last_column = pos;

  return pos;
}

Integer CsvOutputService::
addColumn(String name_column, ConstArrayView<Real> elems)
{
  Integer pos = m_values_csv.dim2Size();
  m_values_csv.resize(m_values_csv.dim1Size(), pos+1);

  m_name_columns.add(name_column);
  m_size_columns.add(0);

  addElemsColumn(pos, elems);

  return pos;
}

bool CsvOutputService::
addColumns(StringConstArrayView name_columns)
{
  Integer size = name_columns.size();
  if(size == 0) return true;

  Integer pos = m_values_csv.dim2Size();
  m_values_csv.resize(m_values_csv.dim1Size(), pos+size);

  m_name_columns.addRange(name_columns);
  m_size_columns.addRange(IntegerUniqueArray(size, 0));

  m_last_column = pos;

  return true;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

bool CsvOutputService::
addElemRow(Integer pos, Real elem)
{
  if(pos < 0 || pos >= m_values_csv.dim1Size()) return false;

  ArrayView<Real> view = m_values_csv[pos];
  Integer size_row = m_size_rows[pos];

  if(m_values_csv.dim2Size() < size_row + 1) return false;

  view[size_row] = elem;

  m_last_row = pos;
  m_last_column = size_row;

  m_size_rows[pos]++;
  // Il peut y avoir des élements sur la ligne d'après à la même colonne.
  // Exemple : addElemRow(pos=L01, elem=NEW):
  // aaa|C00|C01|C02
  // L00|123|456|789
  // L01|147|NEW|
  // L02|159|753|852
  // Il y a 753 donc la taille de la colonne reste égale à 3.
  m_size_columns[size_row] = std::max(pos+1, m_size_columns[size_row]);

  return true;
}

bool CsvOutputService::
addElemRow(String name_row, Real elem, bool create_if_not_exist)
{
  std::optional<Integer> pos = m_name_rows.span().findFirst(name_row);

  if(pos)                       return addElemRow(pos.value(), elem);
  else if(create_if_not_exist)  return addElemRow(addRow(name_row), elem);
  else                          return false;
}

bool CsvOutputService::
addElemSameRow(Real elem)
{
  if(m_last_row == -1 || m_last_column == -1) return false;
  return addElemRow(m_last_row, elem);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

bool CsvOutputService::
addElemsRow(Integer pos, ConstArrayView<Real> elems)
{
  if(pos < 0 || pos >= m_values_csv.dim1Size()) return false;

  ArrayView<Real> view = m_values_csv[pos];
  Integer size_row = m_size_rows[pos];
  Integer min_size = (elems.size() <= m_values_csv.dim2Size()-size_row 
                      ? elems.size() 
                      : m_values_csv.dim2Size()-size_row);

  for(Integer i = 0; i < min_size; i++) {
    view[i+size_row] = elems[i];
    m_size_columns[i+size_row] = std::max(pos+1, m_size_columns[i+size_row]);
  }
  m_size_rows[pos] += min_size;

  m_last_row = pos;
  m_last_column = m_size_rows[pos]-1;

  return elems.size() <= m_values_csv.dim2Size()-size_row;
}

bool CsvOutputService::
addElemsRow(String name_row, ConstArrayView<Real> elems, bool create_if_not_exist)
{
  std::optional<Integer> pos = m_name_rows.span().findFirst(name_row);

  if(pos)                       return addElemsRow(pos.value(), elems);
  // Permet d'avoir un return bool (sinon on pourrait simplement faire addRow(name_row, elems)).
  else if(create_if_not_exist)  return addElemsRow(addRow(name_row), elems); 
  else                          return false;
}

bool CsvOutputService::
addElemsSameRow(ConstArrayView<Real> elems)
{
  if(m_last_row == -1 || m_last_column == -1) return false;
  return addElemsRow(m_last_row, elems);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

bool CsvOutputService::
addElemColumn(Integer pos, Real elem)
{
  if(pos < 0 || pos >= m_values_csv.dim2Size()) return false;

  Integer size_column = m_size_columns[pos];

  if(m_values_csv.dim1Size() < size_column + 1) return false;

  m_values_csv[size_column][pos] = elem;

  m_last_column = pos;
  m_last_row = size_column;

  m_size_columns[pos]++;
  m_size_rows[size_column] = std::max(pos+1, m_size_rows[size_column]);

  return true;
}

bool CsvOutputService::
addElemColumn(String name_column, Real elem, bool create_if_not_exist)
{
  std::optional<Integer> pos = m_name_columns.span().findFirst(name_column);

  if(pos)                       return addElemColumn(pos.value(), elem);
  else if(create_if_not_exist)  return addElemColumn(addColumn(name_column), elem);
  else                          return false;
}

bool CsvOutputService::
addElemSameColumn(Real elem)
{
  if(m_last_row == -1 || m_last_column == -1) return false;
  return addElemColumn(m_last_column, elem);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

bool CsvOutputService::
addElemsColumn(Integer pos, ConstArrayView<Real> elems)
{
  if(pos < 0 || pos >= m_values_csv.dim2Size()) return false;

  Integer size_column = m_size_columns[pos];
  Integer min_size = (elems.size() <= m_values_csv.dim1Size()-size_column 
                      ? elems.size() 
                      : m_values_csv.dim1Size()-size_column);

  for(Integer i = 0; i < min_size; i++) {
    m_values_csv[i+size_column][pos] = elems[i];
    m_size_rows[i+size_column] = std::max(pos+1, m_size_rows[i+size_column]);
  }
  m_size_columns[pos] += min_size;

  m_last_column = pos;
  m_last_row = m_size_columns[pos] - 1;

  return elems.size() <= m_values_csv.dim1Size()-size_column;
}

bool CsvOutputService::
addElemsColumn(String name_column, ConstArrayView<Real> elems, bool create_if_not_exist)
{
  std::optional<Integer> pos = m_name_columns.span().findFirst(name_column);

  if(pos)                       return addElemsColumn(pos.value(), elems);
  // Permet d'avoir un return bool (sinon on pourrait simplement faire addColumn(name_column, elems)).
  else if(create_if_not_exist)  return addElemsColumn(addColumn(name_column), elems); 
  else                          return false;
}

bool CsvOutputService::
addElemsSameColumn(ConstArrayView<Real> elems)
{
  if(m_last_row == -1 || m_last_column == -1) return false;
  return addElemsColumn(m_last_column, elems);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

bool CsvOutputService::
editElemUp(Real elem)
{
  if(m_last_row == -1 || m_last_column == -1 || m_last_row - 1 < 0) return false;
  m_last_row--;

  // Pas besoin d'ajuster la taille de la colonne car on est sûr que m_size_columns[m_last_column] >= m_last_row.
  if(m_size_rows[m_last_row] <= m_last_column) m_size_rows[m_last_row] = m_last_column+1;

  m_values_csv[m_last_row][m_last_column] = elem;
  return true;
}

bool CsvOutputService::
editElemDown(Real elem)
{
  if(m_last_row == -1 || m_last_column == -1 || m_last_row + 1 >= m_values_csv.dim1Size()) return false;
  m_last_row++;
  
  if(m_size_rows[m_last_row] <= m_last_column) m_size_rows[m_last_row] = m_last_column+1;
  if(m_size_columns[m_last_column] <= m_last_row) m_size_columns[m_last_column] = m_last_row+1;

  m_values_csv[m_last_row][m_last_column] = elem;
  return true;
}

bool CsvOutputService::
editElemLeft(Real elem)
{
  if(m_last_row == -1 || m_last_column == -1 || m_last_column - 1 < 0) return false;
  m_last_column--;

  // Pas besoin d'ajuster la taille de la ligne car on est sûr que m_size_rows[m_last_row] >= m_last_column.
  if(m_size_columns[m_last_column] <= m_last_row) m_size_columns[m_last_column] = m_last_row+1;

  m_values_csv[m_last_row][m_last_column] = elem;
  return true;
}

bool CsvOutputService::
editElemRight(Real elem)
{
  if(m_last_row == -1 || m_last_column == -1 || m_last_column + 1 >= m_values_csv.dim2Size()) return false;
  m_last_column++;

  if(m_size_rows[m_last_row] <= m_last_column) m_size_rows[m_last_row] = m_last_column+1;
  if(m_size_columns[m_last_column] <= m_last_row) m_size_columns[m_last_column] = m_last_row+1;

  m_values_csv[m_last_row][m_last_column] = elem;
  return true;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

bool CsvOutputService::
editElem(Integer pos_x, Integer pos_y, Real elem)
{
  if(pos_x < 0 || pos_x >= m_values_csv.dim2Size() 
  || pos_y < 0 || pos_y >= m_values_csv.dim1Size()) 
    return false;


  if(m_size_columns[pos_x] <= pos_y) m_size_columns[pos_x] = pos_y+1;
  if(m_size_rows[pos_y] <= pos_x) m_size_rows[pos_y] = pos_x+1;


  m_values_csv[pos_y][pos_x] = elem;

  m_last_row = pos_y;
  m_last_column = pos_x;

  return true;
}

bool CsvOutputService::
editElem(String name_column, String name_row, Real elem)
{
  std::optional<Integer> pos_x = m_name_columns.span().findFirst(name_column);
  std::optional<Integer> pos_y = m_name_rows.span().findFirst(name_row);

  if(pos_x && pos_y) return editElem(pos_x.value(), pos_y.value(), elem);
  return false;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real CsvOutputService::
getElem(Integer pos_x, Integer pos_y)
{
  if(pos_x < 0 || pos_x >= m_values_csv.dim2Size() 
  || pos_y < 0 || pos_y >= m_values_csv.dim1Size()) 
    return 0;

  return m_values_csv[pos_y][pos_x];
}

Real CsvOutputService::
getElem(String name_column, String name_row)
{
  std::optional<Integer> pos_x = m_name_columns.span().findFirst(name_column);
  std::optional<Integer> pos_y = m_name_rows.span().findFirst(name_row);

  if(pos_x && pos_y) return getElem(pos_x.value(), pos_y.value());
  return 0;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

RealUniqueArray CsvOutputService::
getRow(Integer pos)
{
  Integer size = getSizeRow(pos);
  RealUniqueArray copie(size);
  for(Integer i = 0; i < size; i++){
    copie[i] = m_values_csv[pos][i];
  }
  return copie;
}

RealUniqueArray CsvOutputService::
getRow(String name_row)
{
  std::optional<Integer> pos_y = m_name_rows.span().findFirst(name_row);
  if(pos_y) return getRow(pos_y.value());
  return RealUniqueArray(0);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

RealUniqueArray CsvOutputService::
getColumn(Integer pos)
{
  Integer size = getSizeColumn(pos);

  RealUniqueArray copie(size);
  for(Integer i = 0; i < size; i++){
    copie[i] = m_values_csv[i][pos];
  }
  return copie;
}

RealUniqueArray CsvOutputService::
getColumn(String name_column)
{
  std::optional<Integer> pos_x = m_name_columns.span().findFirst(name_column);
  if(pos_x) return getColumn(pos_x.value());
  return RealUniqueArray(0);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Integer CsvOutputService::
getSizeRow(Integer pos)
{
  if(pos < 0 || pos >= m_values_csv.dim1Size()) return 0;
  return m_size_rows[pos];
}

Integer CsvOutputService::
getSizeRow(String name_row)
{
  std::optional<Integer> pos_y = m_name_rows.span().findFirst(name_row);
  if(pos_y) return getSizeRow(pos_y.value());
  return 0;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Integer CsvOutputService::
getSizeColumn(Integer pos)
{
  if(pos < 0 || pos >= m_values_csv.dim2Size()) return 0;
  return m_size_columns[pos];
}

Integer CsvOutputService::
getSizeColumn(String name_column)
{
  std::optional<Integer> pos_x = m_name_columns.span().findFirst(name_column);
  if(pos_x) return getSizeColumn(pos_x.value());
  return 0;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Integer CsvOutputService::
getPosRow(String name_row)
{
  std::optional<Integer> pos_y = m_name_rows.span().findFirst(name_row);
  if(pos_y) return pos_y.value();
  return -1;
}

Integer CsvOutputService::
getPosColumn(String name_column)
{
  std::optional<Integer> pos_x = m_name_columns.span().findFirst(name_column);
  if(pos_x) return pos_x.value();
  return -1;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Integer CsvOutputService::
getNumRows()
{
  return m_values_csv.dim1Size();
}

Integer CsvOutputService::
getNumColumns()
{
  return m_values_csv.dim2Size();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Integer CsvOutputService::
addAverageColumn(String name_column)
{
  Integer pos = addColumn(name_column);
  for(Integer i = 0; i < m_values_csv.dim1Size(); i++) {
    Real avg = 0.0;
    ConstArrayView<Real> view = m_values_csv[i];
    for(Integer j = 0; j < view.size()-1; j++) {
      avg += view[j];
    }
    avg /= view.size()-1;
    addElemColumn(pos, avg);
  }
  return pos;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsvOutputService::
print(bool only_P0)
{
  if(only_P0 && mesh()->parallelMng()->commRank() != 0) return;
  pinfo() << "P" << mesh()->parallelMng()->commRank() << " - Ecriture du .csv dans la sortie standard :";
  _print(std::cout);
  pinfo() << "P" << mesh()->parallelMng()->commRank() << " - Fin écriture .csv";
}

void CsvOutputService::
print(Integer only_proc)
{
  if(mesh()->parallelMng()->commRank() != only_proc) return;
  pinfo() << "P" << mesh()->parallelMng()->commRank() << " - Ecriture du .csv dans la sortie standard :";
  _print(std::cout);
  pinfo() << "P" << mesh()->parallelMng()->commRank() << " - Fin écriture .csv";
}

bool CsvOutputService::
writeFile(bool only_P0)
{
  String file_name = _computeFinal();
  bool all_only_P0 = m_path_only_P0 && m_name_tab_only_P0;
  // Le seul cas où tout le monde écrit est si only_P0 == false et all_only_P0 == false.
  if((only_P0 || all_only_P0) && mesh()->parallelMng()->commRank() != 0) return true;


  Directory dir(m_path);
  bool sf = false;
  if(mesh()->parallelMng()->commRank() == 0) {
    sf = dir.createDirectory();
  }
  if(mesh()->parallelMng()->commSize() != 1) {
    sf = mesh()->parallelMng()->reduce(Parallel::ReduceMax, sf?1:0);
  }
  if(sf) return false;


  std::ofstream ofile(dir.file(file_name).localstr());
  if(ofile.fail()) return false;

  _print(ofile);

  ofile.close();
  return true;
}

bool CsvOutputService::
writeFile(String path, bool only_P0)
{
  m_path = path;
  m_path_computed = false;
  return writeFile(only_P0);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

String CsvOutputService::
_computeFinal()
{
  if(!m_path_computed) {
    m_path = _computeAt(m_path, m_path_only_P0);
    m_path_computed = true;
  }
  if(!m_name_tab_computed) {
    m_name_tab = _computeAt(m_name_tab, m_name_tab_only_P0);
    m_name_tab_computed = true;
  }
  return m_name_tab + ".csv";
}

String CsvOutputService::
_computeAt(String name, bool& only_P0)
{
  // On découpe la string là où se trouve les @.
  if(name.startsWith("@")){
    name = "@" + name;
  }
  StringUniqueArray string_splited;
  name.split(string_splited, '@');

  // On traite les mots entre les "@".
  if (string_splited.size() > 1)
  {
    // On recherche "proc_id" dans le tableau (donc @proc_id@ dans le nom).
    std::optional<Integer> proc_id = string_splited.span().findFirst("proc_id");
    // On remplace "@proc_id@" par l'id du proc.
    if (proc_id)
    {
      string_splited[proc_id.value()] = String::fromNumber(mesh()->parallelMng()->commRank());
      only_P0 = false;
    }
    // Il n'y a que P0 qui write.
    else
    {
      only_P0 = true;
    }

    // On recherche "num_procs" dans le tableau (donc @num_procs@ dans le nom).
    std::optional<Integer> num_procs = string_splited.span().findFirst("num_procs");
    // On remplace "@num_procs@" par l'id du proc.
    if (num_procs)
    {
      string_splited[num_procs.value()] = String::fromNumber(mesh()->parallelMng()->commSize());
    }
  }

  // On recombine la chaine.
  String combined = "";
  for (String str : string_splited)
  {
    if(str == "@") continue;
    combined = combined + str;
  }
  return combined;
}

void CsvOutputService::
_print(std::ostream& stream)
{
  stream << std::fixed << m_name_tab << m_separator;

  for(Integer j = 0; j < m_name_columns.size(); j++) {
    stream << m_name_columns[j] << m_separator;
  }
  stream << std::endl;

  for(Integer i = 0; i < m_values_csv.dim1Size(); i++) {
    if(m_name_rows.size() > i) stream << m_name_rows[i];
    stream << m_separator;
    ConstArrayView<Real> view = m_values_csv[i];
    for(Integer j = 0; j < m_values_csv.dim2Size(); j++) {
      stream << view[j] << m_separator;
    }
    stream << std::endl;
  }
}

