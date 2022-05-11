<?xml version="1.0"?>
<case codename="Quicksilver" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>TestMedium</title>
    <timeloop>QAMALoop</timeloop>
  </arcane>

  <meshes>
    <mesh>
      <generator name="Cartesian3D" >
        <face-numbering-version>1</face-numbering-version>

        <nb-part-x>8</nb-part-x> 
        <nb-part-y>8</nb-part-y>
        <nb-part-z>4</nb-part-z>

        <origin>0.0 0.0 0.0</origin>

        <x>
          <length>800.0</length>
          <n>80</n>
        </x>

        <y>
          <length>800.0</length>
          <n>80</n>
        </y>

        <z>
          <length>400.0</length>
          <n>40</n>
        </z>

      </generator>
    </mesh>
  </meshes>

  <q-s>
    <dt>2e-09</dt>
    <boundaryCondition>reflect</boundaryCondition>
    <nSteps>10</nSteps>
    <eMax>20</eMax>
    <eMin>1e-09</eMin>
    <nGroups>230</nGroups>
    <lx>800.0</lx>
    <ly>800.0</ly>
    <lz>400.0</lz>
  </q-s>

  <sampling-m-c>
    <nParticles>25600000</nParticles>
    <lowWeightCutoff>0.001</lowWeightCutoff>
    <fMax>0.1</fMax>
    <seed>1029384756</seed>
  </sampling-m-c>

  <tracking-m-c>
    <particle-exchanger name="BasicParticleExchanger">
      <max-nb-message-without-reduce>-1</max-nb-message-without-reduce>
    </particle-exchanger>
    <geometry>
      <material>sourceMaterial</material>
      <shape>brick</shape>
      <xMax>800</xMax>
      <xMin>0</xMin>
      <yMax>800</yMax>
      <yMin>0</yMin>
      <zMax>400</zMax>
      <zMin>0</zMin>
    </geometry>

    <material>
      <name>sourceMaterial</name>
      <mass>12.011</mass>
      <nIsotopes>20</nIsotopes>
      <nReactions>9</nReactions>
      <sourceRate>1e+10</sourceRate>
      <totalCrossSection>1.5</totalCrossSection>
      <absorptionCrossSection>flat</absorptionCrossSection>
      <fissionCrossSection>flat</fissionCrossSection>
      <scatteringCrossSection>flat</scatteringCrossSection>
      <absorptionCrossSectionRatio>0.04</absorptionCrossSectionRatio>
      <fissionCrossSectionRatio>0.05</fissionCrossSectionRatio>
      <scatteringCrossSectionRatio>1</scatteringCrossSectionRatio>
    </material>

    <cross_section>
      <name>flat</name>
      <A>0</A>
      <B>0</B>
      <C>0</C>
      <D>0</D>
      <E>1</E>
      <nuBar>1.6</nuBar>
    </cross_section>

  </tracking-m-c>

</case>
<!-- 

 -->