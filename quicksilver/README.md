# Quicksilver Arcane MiniApp (QAMA)
## Monte-Carlo particle transport problem using Arcane framework

### How to use it?

Build and Install Arcane (version 3.7.3+).

And:

```sh
ARCANE_INSTALL_PATH= # your install path.

QS_BUILD_TYPE=Release # or Debug
QS_SOURCE_DIR=arcane-benchs/quicksilver # your src path
QS_BUILD_DIR=build_qama # build dir

QS_EXE=${QS_BUILD_DIR}/src/Quicksilver
QS_ARC=${QS_SOURCE_DIR}/data/tests/ExampleFull.arc
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

In .arc, you have:
- For CSV:
  - 'csvDir'  : path for results
  - 'csvName' : name of csv files

In these names, you can put '@proc_id@' (replaced by rank of process) and/or '@n_proc@' (replaced by number of processes).

You can change 'loadBalancingMat' (possibles values: 'true' / 'false') and 'loadBalancingLoop' (possibles values: integer >= 0).

All arc input files examples:

```sh
# (for unbalanced tests, see data/Unbal/README.md)

QS_ARC=${QS_SOURCE_DIR}/data/tests/ExampleFull.arc
QS_ARC=${QS_SOURCE_DIR}/data/tests/ExampleDefault.arc
QS_ARC=${QS_SOURCE_DIR}/data/tests/TestSmall.arc
QS_ARC=${QS_SOURCE_DIR}/data/tests/TestMedium.arc
QS_ARC=${QS_SOURCE_DIR}/data/tests/TestLarge.arc
QS_ARC=${QS_SOURCE_DIR}/data/tests/TestHuge.arc

QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n128.arc
QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n256.arc
QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n512.arc
QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n1024.arc
QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n2048.arc
QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n4096.arc
QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n8192.arc
QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n131072.arc

QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n128.arc
QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n256.arc
QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n512.arc
QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n1024.arc
QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n2048.arc
QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n4096.arc
QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n8192.arc
QS_ARC=${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n131072.arc

QS_ARC=${QS_SOURCE_DIR}/data/CTS2/CTS2_n128.arc
QS_ARC=${QS_SOURCE_DIR}/data/CTS2/CTS2_n256.arc
QS_ARC=${QS_SOURCE_DIR}/data/CTS2/CTS2_n512.arc
QS_ARC=${QS_SOURCE_DIR}/data/CTS2/CTS2_n1024.arc
QS_ARC=${QS_SOURCE_DIR}/data/CTS2/CTS2_n2048.arc
QS_ARC=${QS_SOURCE_DIR}/data/CTS2/CTS2_n4096.arc
QS_ARC=${QS_SOURCE_DIR}/data/CTS2/CTS2_n8192.arc
QS_ARC=${QS_SOURCE_DIR}/data/CTS2/CTS2_n131072.arc

QS_ARC=${QS_SOURCE_DIR}/data/qs_original/allAbsorb.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/allEscape.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/scatteringOnly.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/homogeneousProblem.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/no.collisions.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/noFission.arc
QS_ARC=${QS_SOURCE_DIR}/data/qs_original/NonFlatXC.arc
```

Original Quicksilver is available here: https://github.com/LLNL/Quicksilver