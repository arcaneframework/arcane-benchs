cmake_minimum_required(VERSION 3.16)
project(Quicksilver LANGUAGES CXX)

find_package(Arcane REQUIRED)

add_executable(Quicksilver main.cc QSModule.cc QS_axl.h InitMCModule.cc InitMC_axl.h
CollisionEvent.cc
CoralBenchmark.cc
CycleTracking.cc
DecompositionObject.cc
DirectionCosine.cc
EnergySpectrum.cc
GlobalFccGrid.cc
GridAssignmentObject.cc
InputBlock.cc
MCT.cc
MC_Adjacent_Facet.cc
MC_Base_Particle.cc
MC_Domain.cc
MC_Facet_Crossing_Event.cc
MC_Fast_Timer.cc
MC_Load_Particle.cc
MC_Location.cc
MC_Particle_Buffer.cc
MC_RNG_State.cc
MC_Segment_Outcome.cc
MC_SourceNow.cc
MacroscopicCrossSection.cc
MeshPartition.cc
MonteCarlo.cc
MpiCommObject.cc
NuclearData.cc
Parameters.cc
ParticleVault.cc
ParticleVaultContainer.cc
PopulationControl.cc
SendQueue.cc
SharedMemoryCommObject.cc
Tallies.cc
cmdLineParser.cc
cudaFunctions.cc
initMC.cc
main.cc
parseUtils.cc
utils.cc
utilsMpi.cc 
NuclearDataArc.cc )

arcane_generate_axl(QS)
arcane_generate_axl(InitMC)
arcane_add_arcane_libraries_to_target(Quicksilver)
target_include_directories(Quicksilver PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Quicksilver.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
