<?xml version="1.0"?>
<case codename="QAMA" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>allAbsorb</title>
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
          <length>10.0</length>
          <n>10</n>
        </x>

        <y>
          <length>10.0</length>
          <n>10</n>
        </y>

        <z>
          <length>10.0</length>
          <n>10</n>
        </z>

      </generator>
    </mesh>
  </meshes>

  <q-s>
    <dt>1e6</dt>
    <boundaryCondition>reflect</boundaryCondition>
    <nSteps>20</nSteps>
    <eMax>1</eMax>
    <eMin>0.99999</eMin>
    <nGroups>1</nGroups>
    <csvDir>./allAbsorb/</csvDir>
    <csvName>Results_P@proc_id@</csvName>
  </q-s>

  <sampling-m-c>
    <fMax>0.1</fMax>
    <seed>1029384761</seed>
    <nParticles>9999</nParticles>
  </sampling-m-c>

  <tracking-m-c>
    <geometry>
      <material>boxMaterial</material>
      <shape>brick</shape>
      <xMax>10</xMax>
      <xMin>0</xMin>
      <yMax>10</yMax>
      <yMin>0</yMin>
      <zMax>10</zMax>
      <zMin>0</zMin>
    </geometry>

    <geometry>
      <material>sourceMaterial</material>
      <shape>brick</shape>
      <xMax>1</xMax>
      <xMin>0</xMin>
      <yMax>1</yMax>
      <yMin>0</yMin>
      <zMax>1</zMax>
      <zMin>0</zMin>
    </geometry>

    <material>
      <name>boxMaterial</name>
      <nIsotopes>10</nIsotopes>
      <nReactions>9</nReactions>
      <sourceRate>0</sourceRate>
      <totalCrossSection>1e10</totalCrossSection>
      <absorptionCrossSection>flat</absorptionCrossSection>
      <fissionCrossSection>flat</fissionCrossSection>
      <scatteringCrossSection>flat</scatteringCrossSection>
      <absorptionCrossSectionRatio>1</absorptionCrossSectionRatio>
      <fissionCrossSectionRatio>0</fissionCrossSectionRatio>
      <scatteringCrossSectionRatio>0</scatteringCrossSectionRatio>
    </material>

    <material>
      <name>sourceMaterial</name>
      <nIsotopes>10</nIsotopes>
      <nReactions>9</nReactions>
      <sourceRate>1e-2</sourceRate>
      <totalCrossSection>1e10</totalCrossSection>
      <absorptionCrossSection>flat</absorptionCrossSection>
      <fissionCrossSection>flat</fissionCrossSection>
      <scatteringCrossSection>flat</scatteringCrossSection>
      <absorptionCrossSectionRatio>1</absorptionCrossSectionRatio>
      <fissionCrossSectionRatio>0</fissionCrossSectionRatio>
      <scatteringCrossSectionRatio>0</scatteringCrossSectionRatio>
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
  </tracking-m-c>

</case>
