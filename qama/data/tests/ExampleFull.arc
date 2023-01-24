<?xml version="1.0"?>
<case codename="QAMA" xml:lang="en" codeversion="1.0">
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
          <length>100.0</length>
          <n>4</n>
        </x>

        <y>
          <length>100.0</length>
          <n>4</n>
        </y>

        <z>
          <length>100.0</length>
          <n>4</n>
        </z>

      </generator>
    </mesh>
  </meshes>

  <q-s>
    <!-- <dt>2e-10</dt> -->
    <dt>2e-09</dt>
    <boundaryCondition>reflect</boundaryCondition>
    <nSteps>10</nSteps>
    <eMax>20</eMax>
    <eMin>1e-09</eMin>
    <nGroups>230</nGroups>
    <csvDir>./ExampleFull/</csvDir>
    <csvName>Results_P@proc_id@</csvName>
    <loadBalancingMat>true</loadBalancingMat>
    <loadBalancingLoop>0</loadBalancingLoop>
    <csvReferenceDir>default</csvReferenceDir>
    <csvOverwriteReference>false</csvOverwriteReference>
  </q-s>

  <sampling-m-c>
    <nParticles>100000</nParticles>
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
      <!-- <sourceRate>5</sourceRate> -->
      <sourceRate>1e+10</sourceRate>
      <totalCrossSection>1.5</totalCrossSection>
      <absorptionCrossSection>flat</absorptionCrossSection>
      <fissionCrossSection>flat</fissionCrossSection>
      <scatteringCrossSection>flat</scatteringCrossSection>
      <!-- <absorptionCrossSectionRatio>1</absorptionCrossSectionRatio> -->
      <!-- <fissionCrossSectionRatio>0</fissionCrossSectionRatio> -->
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
*I-QS             Number of particles sourced in                              (m_source): 9984
*I-QS             Number of particles Russian Rouletted in population control     (m_rr): 0
*I-QS             Number of particles split in population control              (m_split): 89999
*I-QS             Number of particles absorbed                                (m_absorb): 20312
*I-QS             Number of scatters                                         (m_scatter): 509005
*I-QS             Number of fission events                                   (m_fission): 25438
*I-QS             Number of particles created by collisions                  (m_produce): 40655
*I-QS             Number of collisions                                     (m_collision): 554755
*I-QS             Number of particles that escape                             (m_escape): 0
*I-QS             Number of particles that enter census                       (m_census): 94888
*I-QS             Number of segements                                   (m_num_segments): 672040
*I-QS             Number of particles at end of cycle                            (m_end): 94888
*I-QS             Particles contribution to the scalar flux      (sum_scalar_flux_tally): 73772081.4358831

*I-QS         End iteration #2
*I-QS           Informations:
*I-QS             Number of particles at beginning of cycle                    (m_start): 94888
*I-QS             Number of particles sourced in                              (m_source): 9984
*I-QS             Number of particles Russian Rouletted in population control     (m_rr): 4831
*I-QS             Number of particles split in population control              (m_split): 0
*I-QS             Number of particles absorbed                                (m_absorb): 36110
*I-QS             Number of scatters                                         (m_scatter): 889113
*I-QS             Number of fission events                                   (m_fission): 43996
*I-QS             Number of particles created by collisions                  (m_produce): 70307
*I-QS             Number of collisions                                     (m_collision): 969219
*I-QS             Number of particles that escape                             (m_escape): 0
*I-QS             Number of particles that enter census                       (m_census): 90242
*I-QS             Number of segements                                   (m_num_segments): 1098441
*I-QS             Number of particles at end of cycle                            (m_end): 90242
*I-QS             Particles contribution to the scalar flux      (sum_scalar_flux_tally): 200813514.61978
  -->