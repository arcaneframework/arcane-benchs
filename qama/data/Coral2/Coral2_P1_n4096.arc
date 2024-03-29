<?xml version="1.0"?>
<case codename="QAMA" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Coral2_P1_n4096</title>
    <timeloop>QAMALoop</timeloop>
  </arcane>

  <meshes>
    <mesh>
      <generator name="Cartesian3D" >
        <face-numbering-version>1</face-numbering-version>

        <nb-part-x>16</nb-part-x> 
        <nb-part-y>16</nb-part-y>
        <nb-part-z>16</nb-part-z>

        <origin>0.0 0.0 0.0</origin>

        <x>
          <length>64.0</length>
          <n>64</n>
        </x>

        <y>
          <length>64.0</length>
          <n>64</n>
        </y>

        <z>
          <length>32.0</length>
          <n>32</n>
        </z>

      </generator>
    </mesh>
  </meshes>

  <q-s>
    <dt>2e-09</dt>
    <boundaryCondition>reflect</boundaryCondition>
    <nSteps>100</nSteps>
    <eMax>20</eMax>
    <eMin>1e-09</eMin>
    <nGroups>230</nGroups>
    <csvDir>./Coral2_P1_n4096/</csvDir>
    <csvName>Results_P@proc_id@</csvName>
  </q-s>

  <sampling-m-c>
    <nParticles>5242880</nParticles>
    <fMax>0.1</fMax>
    <seed>1029384756</seed>
    <lowWeightCutoff>0.001</lowWeightCutoff>
  </sampling-m-c>

  <tracking-m-c>
    
    <geometry>
      <material>sourceMaterial</material>
      <shape>brick</shape>
      <xMax>10000</xMax>
      <xMin>0</xMin>
      <yMax>10000</yMax>
      <yMin>0</yMin>
      <zMax>10000</zMax>
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
