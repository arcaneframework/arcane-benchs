<?xml version="1.0"?>
<case codename="QAMA" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>NonFlatXC</title>
    <timeloop>QAMALoop</timeloop>
  </arcane>

  <meshes>
    <mesh>
      <generator name="Cartesian3D" >
        <face-numbering-version>1</face-numbering-version>

        <nb-part-x>2</nb-part-x> 
        <nb-part-y>2</nb-part-y>
        <nb-part-z>1</nb-part-z>

        <origin>0.0 0.0 0.0</origin>

        <x>
          <length>100.0</length>
          <n>10</n>
        </x>

        <y>
          <length>100.0</length>
          <n>10</n>
        </y>

        <z>
          <length>100.0</length>
          <n>10</n>
        </z>

      </generator>
    </mesh>
  </meshes>

  <q-s>
    <dt>1e-08</dt>
    <boundaryCondition>reflect</boundaryCondition>
    <nSteps>10</nSteps>
    <eMax>20</eMax>
    <eMin>1e-08</eMin>
    <nGroups>230</nGroups>
    <csvDir>./NonFlatXC/</csvDir>
    <csvName>Results_P@proc_id@</csvName>
  </q-s>

  <sampling-m-c>
    <fMax>0.1</fMax>
    <nParticles>1000000</nParticles>
    <seed>1029384756</seed>
    <lowWeightCutoff>0.001</lowWeightCutoff>
  </sampling-m-c>

  <tracking-m-c>
    <geometry>
      <material>sourceMaterial</material>
      <shape>brick</shape>
      <xMax>100</xMax>
      <xMin>0</xMin>
      <yMax>100</yMax>
      <yMin>0</yMin>
      <zMax>100</zMax>
      <zMin>0</zMin>
    </geometry>

    <material>
      <name>sourceMaterial</name>
      <mass>1000.0</mass>
      <nIsotopes>10</nIsotopes>
      <nReactions>9</nReactions>
      <sourceRate>1e+10</sourceRate>
      <totalCrossSection>6</totalCrossSection>
      <absorptionCrossSection>absorb</absorptionCrossSection>
      <fissionCrossSection>fission</fissionCrossSection>
      <scatteringCrossSection>scatter</scatteringCrossSection>
      <absorptionCrossSectionRatio>6e-3</absorptionCrossSectionRatio>
      <fissionCrossSectionRatio>1</fissionCrossSectionRatio>
      <scatteringCrossSectionRatio>5</scatteringCrossSectionRatio>
    </material>

    <material>
      <name>flatMaterial</name>
      <nIsotopes>20</nIsotopes>
      <nReactions>9</nReactions>
      <sourceRate>1e+10</sourceRate>
      <totalCrossSection>1</totalCrossSection>
      <absorptionCrossSection>flat</absorptionCrossSection>
      <fissionCrossSection>flat</fissionCrossSection>
      <scatteringCrossSection>flat</scatteringCrossSection>
      <absorptionCrossSectionRatio>1</absorptionCrossSectionRatio>
      <fissionCrossSectionRatio>1</fissionCrossSectionRatio>
      <scatteringCrossSectionRatio>1</scatteringCrossSectionRatio>
    </material>

    <cross_section>
      <name>flat</name>
      <A>0</A>
      <B>0</B>
      <C>0</C>
      <D>0</D>
      <E>1</E>
      <nuBar>2.4</nuBar>
    </cross_section>

    <cross_section>
      <name>absorb</name>
      <A>0</A>
      <B>0</B>
      <C>0</C>
      <D>-0.5243</D>
      <E>-2.22</E>
    </cross_section>

    <cross_section>
      <name>fission</name>
      <A>0</A>
      <B>0</B>
      <C>0</C>
      <D>-0.342</D>
      <E>0</E>
      <nuBar>2.4</nuBar>
    </cross_section>

    <cross_section>
      <name>scatter</name>
      <A>0</A>
      <B>0</B>
      <C>0</C>
      <D>0</D>
      <E>0.7</E>
    </cross_section>
  </tracking-m-c>

</case>
