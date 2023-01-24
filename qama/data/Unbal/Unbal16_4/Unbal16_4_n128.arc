<?xml version="1.0"?>
<case codename="QAMA" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Unbal16_4_n128</title>
    <timeloop>QAMALoop</timeloop>
  </arcane>

  <meshes>
    <mesh>
      <generator name="Cartesian3D" >
        <face-numbering-version>1</face-numbering-version>

        <nb-part-x>16</nb-part-x> 
        <nb-part-y>4</nb-part-y>
        <nb-part-z>2</nb-part-z>

        <origin>0.0 0.0 0.0</origin>

        <x>
          <length>128.0</length>
          <n>16</n>
        </x>

        <y>
          <length>256.0</length>
          <n>8</n>
        </y>

        <z>
          <length>256.0</length>
          <n>4</n>
        </z>

      </generator>
    </mesh>
  </meshes>

  <q-s>
    <dt>1e-08</dt>
    <boundaryCondition>octant</boundaryCondition>
    <nSteps>15</nSteps>
    <eMax>20</eMax>
    <eMin>1e-09</eMin>
    <nGroups>230</nGroups>
    <csvDir>./Unbal16_4_n128/</csvDir>
    <csvName>Results_P@proc_id@</csvName>
    <loadBalancingMat>true</loadBalancingMat>
    <loadBalancingLoop>0</loadBalancingLoop>
  </q-s>

  <sampling-m-c>
    <nParticles>512000</nParticles>
    <lowWeightCutoff>0.001</lowWeightCutoff>
  </sampling-m-c>

  <tracking-m-c>

    <particle-exchanger name="BasicParticleExchanger">
      <max-nb-message-without-reduce>-1</max-nb-message-without-reduce>
    </particle-exchanger>

    <!-- 
        |S|H|H|S|S|S|S|H|H|H|H|S|S|S|S|S|
        |S|H|H|S|S|S|S|H|H|H|H|S|S|S|S|S|
        |S|H|H|S|S|S|S|H|H|H|H|S|S|S|S|S|
        |S|H|H|S|S|S|S|H|H|H|H|S|S|S|S|S|
    z
    ^
    |__>x
    -->

    <geometry>
      <material>softComputing</material>
      <shape>brick</shape>

      <xMin>0</xMin>
      <xMax>4</xMax>

      <yMin>0</yMin>
      <yMax>512</yMax>

      <zMin>0</zMin>
      <zMax>512</zMax>
    </geometry>

    <geometry>
      <material>hardComputing</material>
      <shape>brick</shape>

      <xMin>4</xMin>
      <xMax>12</xMax>

      <yMin>0</yMin>
      <yMax>512</yMax>

      <zMin>0</zMin>
      <zMax>512</zMax>
    </geometry>

    <geometry>
      <material>softComputing</material>
      <shape>brick</shape>

      <xMin>12</xMin>
      <xMax>28</xMax>

      <yMin>0</yMin>
      <yMax>512</yMax>

      <zMin>0</zMin>
      <zMax>512</zMax>
    </geometry>

    <geometry>
      <material>hardComputing</material>
      <shape>brick</shape>

      <xMin>28</xMin>
      <xMax>44</xMax>

      <yMin>0</yMin>
      <yMax>512</yMax>

      <zMin>0</zMin>
      <zMax>512</zMax>
    </geometry>

    <geometry>
      <material>softComputing</material>
      <shape>brick</shape>

      <xMin>44</xMin>
      <xMax>64</xMax>

      <yMin>0</yMin>
      <yMax>512</yMax>

      <zMin>0</zMin>
      <zMax>512</zMax>
    </geometry>


    <material>
      <name>hardComputing</name>
      <mass>2</mass>

      <nIsotopes>20</nIsotopes>
      <nReactions>9</nReactions>
      <sourceRate>5e+10</sourceRate>
      <totalCrossSection>1.5</totalCrossSection>

      <absorptionCrossSection>flat0</absorptionCrossSection>
      <fissionCrossSection>flat0</fissionCrossSection>
      <scatteringCrossSection>flat0</scatteringCrossSection>

      <absorptionCrossSectionRatio>0.8</absorptionCrossSectionRatio>
      <fissionCrossSectionRatio>1.0</fissionCrossSectionRatio>
      <scatteringCrossSectionRatio>0.5</scatteringCrossSectionRatio>
    </material>

    <cross_section>
      <name>flat0</name>
      <A>0</A>
      <B>0</B>
      <C>0</C>
      <D>0</D>
      <E>1</E>
      <nuBar>2.0</nuBar>
    </cross_section>


    <material>
      <name>softComputing</name>
      <mass>1</mass>

      <nIsotopes>20</nIsotopes>
      <nReactions>9</nReactions>
      <sourceRate>1e+10</sourceRate>
      <totalCrossSection>1.5</totalCrossSection>

      <absorptionCrossSection>flat1</absorptionCrossSection>
      <fissionCrossSection>flat1</fissionCrossSection>
      <scatteringCrossSection>flat1</scatteringCrossSection>
      
      <absorptionCrossSectionRatio>1</absorptionCrossSectionRatio>
      <fissionCrossSectionRatio>0.05</fissionCrossSectionRatio>
      <scatteringCrossSectionRatio>0.05</scatteringCrossSectionRatio>
    </material>

    <cross_section>
      <name>flat1</name>
      <A>0</A>
      <B>0</B>
      <C>0</C>
      <D>0</D>
      <E>1</E>
      <nuBar>0.5</nuBar>
    </cross_section>

  </tracking-m-c>

</case>
