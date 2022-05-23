<?xml version="1.0"?>
<case codename="Quicksilver" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Coral2_P2_n1024</title>
    <timeloop>QAMALoop</timeloop>
  </arcane>

  <meshes>
    <mesh>
      <generator name="Cartesian3D" >
        <face-numbering-version>1</face-numbering-version>

        <nb-part-x>16</nb-part-x> 
        <nb-part-y>8</nb-part-y>
        <nb-part-z>8</nb-part-z>

        <origin>0.0 0.0 0.0</origin>

        <x>
          <length>22.0</length>
          <n>22</n>
        </x>

        <y>
          <length>22.0</length>
          <n>22</n>
        </y>

        <z>
          <length>22.0</length>
          <n>22</n>
        </z>

      </generator>
    </mesh>
  </meshes>

  <q-s>
    <dt>1e-08</dt>
    <boundaryCondition>reflect</boundaryCondition>
    <nSteps>100</nSteps>
    <eMax>20</eMax>
    <eMin>1e-08</eMin>
    <nGroups>230</nGroups>
    <csvFile>./Coral2_P2_n1024.csv</csvFile>
  </q-s>

  <sampling-m-c>
    <nParticles>425920</nParticles>
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
      <mass>1.5</mass>
      <nIsotopes>10</nIsotopes>
      <nReactions>3</nReactions>
      <sourceRate>1e+10</sourceRate>
      <totalCrossSection>16.75</totalCrossSection>
      <absorptionCrossSection>absorb</absorptionCrossSection>
      <fissionCrossSection>fission</fissionCrossSection>
      <scatteringCrossSection>scatter</scatteringCrossSection>
      <absorptionCrossSectionRatio>10</absorptionCrossSectionRatio>
      <fissionCrossSectionRatio>8</fissionCrossSectionRatio>
      <scatteringCrossSectionRatio>82</scatteringCrossSectionRatio>
    </material>

    <cross_section>
      <name>absorb</name>
      <A>0</A>
      <B>0</B>
      <C>0</C>
      <D>-0.2</D>
      <E>2</E>
    </cross_section>

    <cross_section>
      <name>fission</name>
      <A>0</A>
      <B>0</B>
      <C>0</C>
      <D>-0.2</D>
      <E>2</E>
      <nuBar>2</nuBar>
    </cross_section>

    <cross_section>
      <name>scatter</name>
      <A>0</A>
      <B>0</B>
      <C>0</C>
      <D>0</D>
      <E>97</E>
    </cross_section>

  </tracking-m-c>

</case>
