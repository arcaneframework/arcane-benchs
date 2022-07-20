<?xml version="1.0"?>
<case codename="Synchronize" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Test des synchronisations</title>
    <timeloop>SynchronizeLoop</timeloop>
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
</case>
