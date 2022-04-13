<?xml version="1.0"?>
<case codename="Quicksilver" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>QuicksilverLoop</timeloop>
  </arcane>

  <meshes>
    <mesh>
      <generator name="Cartesian3D" >
        <face-numbering-version>1</face-numbering-version>

        <nb-part-x>1</nb-part-x> 
        <nb-part-y>1</nb-part-y>
        <nb-part-z>1</nb-part-z>

        <origin>0.0 0.0 0.0</origin>

        <!-- lx = lenght * n -->
        <x>
          <length>10.0</length>
          <n>3</n>
        </x>

        <!-- ly = lenght * n -->
        <y>
          <length>10.0</length>
          <n>3</n>
        </y>

        <!-- lz = lenght * n -->
        <z>
          <length>10.0</length>
          <n>3</n>
        </z>

      </generator>
      <!-- <initialisation>
        <variable nom="Tally" valeur="0."/>
      </initialisation> -->
    </mesh>
  </meshes>

  <qt>
    <dt>2e-09</dt>
    <fMax>0.1</fMax>
    <!-- <nParticles>1000000</nParticles> -->
    <nSteps>1</nSteps>
    <seed>1029384756</seed>
    <lx>100.0</lx>
    <ly>100.0</ly>
    <lz>100.0</lz>
    <boundaryCondition>reflect</boundaryCondition>
    <eMax>20</eMax>
    <eMin>1e-09</eMin>
    <nGroups>230</nGroups>
    <lowWeightCutoff>0.001</lowWeightCutoff>

    <!-- <material>sourceMaterial</material>
    <shape>brick</shape>
    <xMax>10000</xMax>
    <xMin>0</xMin>
    <yMax>10000</yMax>
    <yMin>0</yMin>
    <zMax>10000</zMax>
    <zMin>0</zMin> -->

  </qt>
</case>
