// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "ArcaneLoops_axl.h"
#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>

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

 private:

  Timer m_timer;
  CellGroup m_cells;

 private:

  void _doCellNodesViaNodes();
  void _doCellNodesDirect();
  void _doNodesCellNodes();
  void _doNodesCellNodesDirect();
  void _doNodesCellNodesDirect2();
  void _doTest(void (ArcaneLoopsModule::*functor)(), const String& name,Integer nb_z = 20000);
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
, m_timer(mbi.subDomain(),"Test",Timer::TimerReal)
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneLoopsModule::
compute()
{
  info() << "Module ArcaneLoops COMPUTE";

  _doTest(&ArcaneLoopsModule::_doCellNodesViaNodes,"CellNodesViaNodes");
  _doTest(&ArcaneLoopsModule::_doCellNodesDirect,"CellNodesDirect");
  _doTest(&ArcaneLoopsModule::_doNodesCellNodes,"NodesCellNodes",4000);
  _doTest(&ArcaneLoopsModule::_doNodesCellNodesDirect,"NodesCellNodesDirect",4000);
  _doTest(&ArcaneLoopsModule::_doNodesCellNodesDirect2,"NodesCellNodesDirect2",4000);

  // Stop code after 10 iterations
  //if (m_global_iteration()>10)
  subDomain()->timeLoopMng()->stopComputeLoop(true);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneLoopsModule::
startInit()
{
  info() << "Module ArcaneLoops INIT";
  m_cells = allCells();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneLoopsModule::
_doCellNodesViaNodes()
{
  static int nb_iter = 0;
  Int64 total = 0;
  ENUMERATE_CELL(icell,m_cells){
    const Cell& cell = *icell;
    Integer nb_node = cell.nbNode();
    for( Integer i=0; i<nb_node; ++i ){
      const Node& n = cell.node(i);
      total += n.uniqueId().asInt64();
    }
  }
  _setTotal(total);
  if (nb_iter==0)
    info() << "TOTAL=" << total;
  ++nb_iter;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneLoopsModule::
_doCellNodesDirect()
{
  static int nb_iter = 0;
  Int64 total = 0;
  ENUMERATE_CELL(icell,m_cells){
    const Cell& cell = *icell;
    for( Node n : cell.nodes() ){
      total += n.uniqueId().asInt64();
    }
  }
  _setTotal(total);
  if (nb_iter==0)
    info() << "TOTAL=" << total;
  ++nb_iter;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneLoopsModule::
_doNodesCellNodes()
{
  static int nb_iter = 0;
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
  _setTotal(total);
  if (nb_iter==0)
    info() << "TOTAL1=" << total;
  ++nb_iter;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneLoopsModule::
_doNodesCellNodesDirect()
{
  static int nb_iter = 0;
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
  _setTotal(total);
  if (nb_iter==0)
    info() << "TOTAL1=" << total;
  ++nb_iter;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneLoopsModule::
_doNodesCellNodesDirect2()
{
  static int nb_iter = 0;
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
  _setTotal(total);
  if (nb_iter==0)
    info() << "TOTAL1=" << total;
  ++nb_iter;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneLoopsModule::
_doTest( void (ArcaneLoopsModule::*functor)(), const String& name, Integer nb_z)
{
  info() << "Begin test name=" << name;
  // Réduit la taille du test en débug pour qu'il ne dure pas trop longtemps
  if (arcaneIsDebug())
    nb_z /= 20;

  {
    Timer::Sentry ts(&m_timer);
    for( Integer k=0; k<nb_z; ++k ){
      (this->*functor)();
    }
  }
  info()  << " name=" << Trace::Width(30) << name << " TIME = " << m_timer.lastActivationTime();
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
