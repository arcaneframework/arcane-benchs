<?xml version="1.0"?>
<case codename="Synchronize" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Test des synchronisations</title>
    <timeloop>SynchronizeLoop</timeloop>
    <modules>
      <module name="ArcaneLoadBalance" active="true" />
    </modules>
  </arcane>

  <meshes>
    <mesh>
      <ghost-layer-builder-version>4</ghost-layer-builder-version>
      <generator name="Cartesian3D" >
        <nb-part-x>2</nb-part-x> 
        <nb-part-y>2</nb-part-y>
        <nb-part-z>2</nb-part-z>
        <origin>1.0 2.0 3.0</origin>
        <x><n>100</n><length>1.0</length></x>
        <y><n>20</n><length>0.3</length></y>
        <z><n>20</n><length>0.3</length></z>
      </generator>
    </mesh>
  </meshes>

  <arcane-load-balance>
    <active>true</active>
    <partitioner name="MeshPartitionerTester">
      <sub-rank-divider>8</sub-rank-divider>
    </partitioner>
    <period>5</period>
    <statistics>true</statistics>
    <max-imbalance>0.01</max-imbalance>
    <min-cpu-time>0</min-cpu-time>
  </arcane-load-balance>

  <synchronize>
    <nb-synchronize-per-iteration>20000</nb-synchronize-per-iteration>
  </synchronize>

</case>
