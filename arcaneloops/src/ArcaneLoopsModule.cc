// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "ArcaneLoops_axl.h"
#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/UnstructuredMeshConnectivity.h>

using namespace Arcane;

/*!
 * \brief Module ArcaneLoops.
 */
class ArcaneLoopsModule
: public ArcaneArcaneLoopsObject
{
 public:
  explicit ArcaneLoopsModule(const ModuleBuildInfo& mbi);

 public:
  /*!
   * \brief Méthode appelée à chaque itération.
   */
  void compute() override;
  /*!
   * \brief Méthode appelée lors de l'initialisation.
   */
  void startInit() override;

  /** Retourne le numéro de version du module */
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

 private:

  VariableCellReal m_var1;
  VariableCellReal m_var2;
  VariableCellReal m_var3;
  VariableCellReal m_var4;
  VariableCellReal m_var5;
  VariableNodeInt64 m_node_unique_ids;
  VariableCellInt64 m_cell_unique_ids;
  Int64 m_full_total = 0;

 private:

  Timer m_timer;
  CellGroup m_cells;
  UnstructuredMeshConnectivityView m_connectivity_view;

 private:

  Int64 _doCellNodesViaNodes();
  Int64 _doCellNodesDirect();
  Int64 _doNodesCellNodes();
  Int64 _doNodesCellNodesDirect();
  Int64 _doNodesCellNodesDirect2();
  Int64 _doNodesCellNodesWithConnectivity();
  Int64 _doCellLocalIdFromGroup2(CellGroup cells);
  Int64 _doCellLocalIdFromItemVector2(CellVectorView cells);
  Int64 _doCellUniqueIdFromGroup2(CellGroup cells);
  Int64 _doCellUniqueIdFromItemVector2(CellVectorView cells);
  Int64 _doCellOperation1(CellGroup cells);

  Int64 _doCellLocalIdFromGroup() { return _doCellLocalIdFromGroup2(m_cells); }
  Int64 _doCellLocalIdFromItemVector() { return _doCellLocalIdFromItemVector2(m_cells.view()); }
  Int64 _doCellUniqueIdFromGroup() { return _doCellUniqueIdFromGroup2(m_cells); }
  Int64 _doCellUniqueIdFromItemVector() { return _doCellUniqueIdFromItemVector2(m_cells.view()); }
  Int64 _doCellOperation1X() { return _doCellOperation1(m_cells); }

  void _doTest(const String& id,Int64 (ArcaneLoopsModule::*functor)(), const String& name,Integer nb_z = 20000);
  void _setTotal(Int64 total);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ArcaneLoopsModule::
ArcaneLoopsModule(const ModuleBuildInfo& mbi) 
: ArcaneArcaneLoopsObject(mbi)
, m_var1(VariableBuildInfo(this,"Var1"))
, m_var2(VariableBuildInfo(this,"Var2"))
, m_var3(VariableBuildInfo(this,"Var3"))
, m_var4(VariableBuildInfo(this,"Var4"))
, m_var5(VariableBuildInfo(this,"Var5"))
, m_node_unique_ids(VariableBuildInfo(this,"NodeUniqueIds"))
, m_cell_unique_ids(VariableBuildInfo(this,"CellUniqueIds"))
, m_timer(mbi.subDomain(),"Test",Timer::TimerReal)
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneLoopsModule::
compute()
{
  info() << "Module ArcaneLoops COMPUTE";
  m_full_total = 0;

  _doTest("001",&ArcaneLoopsModule::_doCellNodesViaNodes,"CellNodesViaNodes");
  _doTest("002",&ArcaneLoopsModule::_doCellNodesDirect,"CellNodesDirect");
  _doTest("003",&ArcaneLoopsModule::_doNodesCellNodes,"NodesCellNodes",4000);
  _doTest("004",&ArcaneLoopsModule::_doNodesCellNodesDirect,"NodesCellNodesDirect",4000);
  _doTest("005",&ArcaneLoopsModule::_doNodesCellNodesDirect2,"NodesCellNodesDirect2",4000);
  _doTest("006",&ArcaneLoopsModule::_doCellLocalIdFromGroup,"CellLocalIdFromGroup",50000);
  _doTest("007",&ArcaneLoopsModule::_doCellLocalIdFromItemVector,"CellLocalIdFromItemVector",50000);
  _doTest("008",&ArcaneLoopsModule::_doCellUniqueIdFromGroup,"CellUniqueIdFromGroup",50000);
  _doTest("009",&ArcaneLoopsModule::_doCellUniqueIdFromItemVector,"CellUniqueIdFromItemVector",50000);
  _doTest("010",&ArcaneLoopsModule::_doCellOperation1X,"CellOperation1",30000);

  subDomain()->timeLoopMng()->stopComputeLoop(true);
  info() << "FULL_TOTAL=" << m_full_total;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneLoopsModule::
startInit()
{
  info() << "Module ArcaneLoops INIT";
  m_cells = allCells();
  m_connectivity_view.setMesh(this->mesh());
  ENUMERATE_(Node,inode,allNodes()){
    m_node_unique_ids[inode] = (*inode).uniqueId();
  }
  ENUMERATE_(Cell,icell,allCells()){
    Cell cell = *icell;
    Int32 lid = cell.localId();
    m_cell_unique_ids[icell] = cell.uniqueId();
    Real x = (Real)(lid+1);
    m_var1[cell] = x*1.2;
    m_var2[cell] = x*1.5;
    m_var3[cell] = -x*1.1;
    m_var4[cell] = x*2.3;
    m_var5[cell] = 0.0;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Int64 ArcaneLoopsModule::
_doCellNodesViaNodes()
{
  Int64 total = 0;
  ENUMERATE_CELL(icell,m_cells){
    const Cell& cell = *icell;
    Integer nb_node = cell.nbNode();
    for( Integer i=0; i<nb_node; ++i ){
      const Node& n = cell.node(i);
      total += n.uniqueId().asInt64();
    }
  }
  return total;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Int64 ArcaneLoopsModule::
_doCellNodesDirect()
{
  Int64 total = 0;
  ENUMERATE_CELL(icell,m_cells){
    const Cell& cell = *icell;
    for( Node n : cell.nodes() ){
      total += n.uniqueId().asInt64();
    }
  }
  return total;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Int64 ArcaneLoopsModule::
_doNodesCellNodes()
{
  Int64 total = 0;
  ENUMERATE_NODE(inode,allNodes()){
    const Node& node = *inode;
    CellEnumerator node_cells = node.cells();
    for( CellEnumerator icell(node_cells); icell(); ++icell ){      
      const Cell& cell = *icell;
      const Integer nb_node = cell.nbNode();
      for( Integer i=0; i<nb_node; ++i ){
        const Node& local_node = cell.node(i);
        if (local_node.uniqueId()==node.uniqueId())
          total += local_node.uniqueId().asInt64();
      }
    }
  }
  return total;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Int64 ArcaneLoopsModule::
_doNodesCellNodesDirect()
{
  Int64 total = 0;
  ENUMERATE_NODE(inode,allNodes()){
    Node node = *inode;
    for( Cell cell : node.cells() ){
      for( Node local_node : cell.nodes() ){
        if (local_node.uniqueId()==node.uniqueId())
          total += local_node.uniqueId().asInt64();
      }
    }
  }
  return total;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Int64 ArcaneLoopsModule::
_doNodesCellNodesDirect2()
{
  Int64 total = 0;
  ENUMERATE_NODE(inode,allNodes()){
    Node node = *inode;
    for( CellEnumerator icell(node.cells()); icell(); ++icell ){      
      Cell cell = *icell;
      const Integer nb_node = cell.nbNode();
      for( NodeEnumerator inode2(cell.nodes()); inode2(); ++inode2 ){      
        Node local_node = *inode2;
        if (local_node.uniqueId()==node.uniqueId())
          total += local_node.uniqueId().asInt64();
      }
    }
  }
  return total;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Int64 ArcaneLoopsModule::
_doNodesCellNodesWithConnectivity()
{
  Int64 total = 0;
  auto c_cell_node = m_connectivity_view.cellNode();
  auto c_node_cell = m_connectivity_view.nodeCell();

  ENUMERATE_NODE(inode,allNodes()){
    NodeLocalId node = *inode;
    for( CellLocalId cell : c_node_cell.cells(node) ){
      for( NodeLocalId local_node : c_cell_node.nodes(cell) ){
        if (m_node_unique_ids[local_node]==m_node_unique_ids[node])
          total += m_node_unique_ids[local_node];
      }
    }
  }
  return total;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Int64 ArcaneLoopsModule::
_doCellLocalIdFromGroup2(CellGroup cells)
{
  Int64 total = 0;
  ENUMERATE_CELL(icell,cells){
    const Cell& cell = *icell;
    total += cell.localId();
  }
  return total;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Int64 ArcaneLoopsModule::
_doCellLocalIdFromItemVector2(CellVectorView cells)
{
  Int64 total = 0;
  for( Cell cell : cells ){
    total += cell.localId();
  }
  return total;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Int64 ArcaneLoopsModule::
_doCellUniqueIdFromGroup2(CellGroup cells)
{
  Int64 total = 0;
  ENUMERATE_CELL(icell,cells){
    const Cell& cell = *icell;
    total += cell.uniqueId().asInt64();
  }
  return total;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Int64 ArcaneLoopsModule::
_doCellUniqueIdFromItemVector2(CellVectorView cells)
{
  Int64 total = 0;
  for( Cell cell : cells ){
    total += cell.uniqueId().asInt64();
  }
  return total;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Int64 ArcaneLoopsModule::
_doCellOperation1(CellGroup cells)
{
  ENUMERATE_(Cell,icell,cells){
    //Cell cell = *icell;
    //m_var5[cell] = m_var2[cell] * m_var3[cell] + m_var4[cell];
    m_var5[icell] = m_var2[icell] * m_var3[icell] + m_var4[icell];
  }
  return 0;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneLoopsModule::
_doTest(const String& id,Int64 (ArcaneLoopsModule::*functor)(), const String& name, Integer nb_z)
{
  info(4) << "Begin test name=" << name;

  // Réduit la taille du test en débug pour qu'il ne dure pas trop longtemps
  if (arcaneIsDebug())
    nb_z /= 20;
  Int32 multiplier = options()->nbLoopMultiplier();
  nb_z *= multiplier;
  nb_z /= 100;

  Int64 full_total = 0;
  {
    Timer::Sentry ts(&m_timer);
    for( Integer k=0; k<nb_z; ++k ){
      full_total += (this->*functor)();
    }
  }
  info()  << "id=" << id << " name=" << Trace::Width(30) << name << " TIME = " << m_timer.lastActivationTime();
  m_full_total += full_total;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneLoopsModule::
_setTotal(Int64 total)
{
  auto x = defaultMesh()->cellFamily()->itemsInternal();
  Cell cell(x[0]);
  m_var1[cell] = (Real)total;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_ARCANELOOPS(ArcaneLoopsModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
