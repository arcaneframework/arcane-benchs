# Quicksilver Arcane MiniApp (QAMA)
## Monte-Carlo particle transport problem using Arcane framework

### How to use it?

Build and Install Arcane.

And:

```sh
ARCANE_INSTALL_PATH= # your install path.

QS_BUILD_TYPE=Release # or Debug
QS_SOURCE_DIR=arcane-benchs/quicksilver # your src path
QS_BUILD_DIR=build_qama # build dir

QS_EXE=${QS_BUILD_DIR}/src/Quicksilver
QS_ARC=${QS_SOURCE_DIR}/data/ExampleFull.arc
cd ${QS_BUILD_DIR}

cmake -S ${QS_SOURCE_DIR} -B ${QS_BUILD_DIR} -DCMAKE_PREFIX_PATH=${ARCANE_INSTALL_PATH} -DCMAKE_BUILD_TYPE=${QS_BUILD_TYPE}
make

# Sequential:
${QS_EXE} ${QS_ARC}

# MPI 4 procs:
mpirun -n 4 ${QS_EXE} ${QS_ARC}

# 2 tasks:
${QS_EXE} -A,T=2 ${QS_ARC}

# MPI 4 procs and 2 tasks/proc:
mpirun -n 4 ${QS_EXE} -A,T=2 ${QS_ARC}
```

Arc input files examples:

```sh
QS_ARC=${QS_SOURCE_DIR}/data/ExampleFull.arc
QS_ARC=${QS_SOURCE_DIR}/data/ExampleDefault.arc
QS_ARC=${QS_SOURCE_DIR}/data/TestSmall.arc
QS_ARC=${QS_SOURCE_DIR}/data/TestMedium.arc
QS_ARC=${QS_SOURCE_DIR}/data/TestLarge.arc
QS_ARC=${QS_SOURCE_DIR}/data/TestHuge.arc

QS_ARC=${QS_SOURCE_DIR}/data/qs_original/Coral2_P1.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/Coral2_P2.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/allAbsorb.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/allEscape.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/scatteringOnly.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/CTS2.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/homogeneousProblem.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/no.collisions.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/noFission.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/NonFlatXC.arc
```

Original Quicksilver is available here: https://github.com/LLNL/Quicksilver