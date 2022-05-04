<?xml version="1.0"?>
<case codename="Quicksilver" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>ExampleFull</title>
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
          <n>4</n>
        </x>

        <y>
          <length>10.0</length>
          <n>4</n>
        </y>

        <z>
          <length>10.0</length>
          <n>4</n>
        </z>

      </generator>
    </mesh>
  </meshes>

  <q-s>
    <dt>2e-09</dt>
    <boundaryCondition>reflect</boundaryCondition>
    <nSteps>2</nSteps>
    <eMax>20</eMax>
    <eMin>1e-09</eMin>
    <nGroups>230</nGroups>
    <lx>100.0</lx>
    <ly>100.0</ly>
    <lz>100.0</lz>
  </q-s>

  <sampling-m-c>
    <nParticles>1000000</nParticles>
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
<!-- 
*I-QS         End iteration #1
*I-QS           Informations:
*I-QS             Number of particles at beginning of cycle                    (m_start): 0
*I-QS             Number of particles sourced in                              (m_source): 99968
*I-QS             Number of particles Russian Rouletted in population control     (m_rr): 0
*I-QS             Number of particles split in population control              (m_split): 900005
*I-QS             Number of particles absorbed                                (m_absorb): 203417
*I-QS             Number of scatters                                         (m_scatter): 5082229
*I-QS             Number of fission events                                   (m_fission): 254032
*I-QS             Number of particles created by collisions                  (m_produce): 406826
*I-QS             Number of collisions                                     (m_collision): 5539678
*I-QS             Number of particles that escape                             (m_escape): 0
*I-QS             Number of particles that enter census                       (m_census): 949350
*I-QS             Number of segements                                   (m_num_segments): 6711690
*I-QS             Number of particles at end of cycle                            (m_end): 949350

*I-QS         End iteration #2
*I-QS           Informations:
*I-QS             Number of particles at beginning of cycle                    (m_start): 949350
*I-QS             Number of particles sourced in                              (m_source): 99968
*I-QS             Number of particles Russian Rouletted in population control     (m_rr): 49354
*I-QS             Number of particles split in population control              (m_split): 0
*I-QS             Number of particles absorbed                                (m_absorb): 356811
*I-QS             Number of scatters                                         (m_scatter): 8909894
*I-QS             Number of fission events                                   (m_fission): 445311
*I-QS             Number of particles created by collisions                  (m_produce): 711828
*I-QS             Number of collisions                                     (m_collision): 9712016
*I-QS             Number of particles that escape                             (m_escape): 0
*I-QS             Number of particles that enter census                       (m_census): 909670
*I-QS             Number of segements                                   (m_num_segments): 11013245
*I-QS             Number of particles at end of cycle                            (m_end): 909670
 -->