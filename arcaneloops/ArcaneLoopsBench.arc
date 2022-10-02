<?xml version="1.0"?>
<case codename="ArcaneLoopsBench" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>ArcaneLoopsBenchLoop</timeloop>
  </arcane>

  <meshes>
    <mesh>
      <generator name="Cartesian2D" >
        <nb-part-x>1</nb-part-x> 
        <nb-part-y>1</nb-part-y>
        <origin>0.0 0.0</origin>
        <x><n>100</n><length>2.0</length></x>
        <y><n>100</n><length>2.0</length></y>
      </generator>
    </mesh>
  </meshes>

  <arcane-loops>
    <nb-loop-multiplier>100</nb-loop-multiplier>
  </arcane-loops>
</case>
