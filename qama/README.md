# Quicksilver Arcane MiniApp (QAMA)
### Monte-Carlo particle transport problem using Arcane framework

## How to use it?

Build and Install Arcane (version 3.7.3+).

And, if you want to compile only QAMA:

```sh
ARCANE_INSTALL_PATH= # your install path.

QS_BUILD_TYPE=Release # or Debug
QS_SOURCE_DIR=arcane-benchs/qama # your src path
QS_BUILD_DIR=build_qama # build dir

QS_EXE=${QS_BUILD_DIR}/src/qama
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

## Output results in CSV file

QAMA uses SimpleCsvOutput to provide an easy method for outputting results
to CSV files.

In the .arc file, you have the following options:
- `csvDir`  : directory for the output CSV files
- `csvName` : name of the output CSV files

In these file names, you can use `@proc_id@` (which will be replaced by the rank of the process) and/or `@n_proc@` (which will be replaced by the number of processes) for customizing the file names.

Example:
```xml
<csvDir>ExampleFull</csvDir>
<csvName>Results_P@proc_id@</csvName>
```
With this example, the CSV files will be located in `./output/csv/ExampleFull/Results_P0.csv`.


## Compare results with CSV reference file

QAMA uses SimpleCsvComparator to provide a method for comparing the current
results with a reference result. Reference files can be created by you,
or you can use the reference files available in this repository.

To use SimpleCsvComparator, you have two options in the .arc file:
```xml
<csvReferenceDir>default</csvReferenceDir>
<csvOverwriteReference>false</csvOverwriteReference>
```
The csvReferenceDir option specifies the directory that contains all the reference files.
SimpleCsvComparator uses the options of SimpleCsvOutput to determine the exact path of the directory.

For example:
```xml
<csvReferenceDir>default</csvReferenceDir>
<csvDir>ExampleFull</csvDir>
<csvName>Results_P@proc_id@</csvName>
```
The path determined by SimpleCsvComparator is :
`./output/csv_refs/ExampleFull/Results_P0.csv`

```xml
<csvReferenceDir>/tmp/my_tmp_refs</csvReferenceDir>
<csvDir>Example123</csvDir>
<csvName>Results456_P@proc_id@</csvName>
```
The path determined by SimpleCsvComparator is:
`tmp/my_tmp_refs/Example123/Results456_P0.csv`

```xml
<csvReferenceDir>/tmp/my_tmp_refs</csvReferenceDir>
<!-- <csvDir></csvDir> -->
<csvName>Results456_P@proc_id@</csvName>
```
The path determined by SimpleCsvComparator is :
`tmp/my_tmp_refs/Results456_P0.csv`

```xml
<csvReferenceDir>/tmp/my_tmp_refs</csvReferenceDir>
<!-- <csvDir></csvDir> -->
<!-- <csvName></csvName> -->
```
The path determined by SimpleCsvComparator is :
`tmp/my_tmp_refs/Table_P0.csv`

```xml
<!-- <csvReferenceDir></csvReferenceDir> -->
```
No comparison will be made.

___

You also have the csvOverwriteReference option. This option is used to write
or overwrite the reference files in the csvDir directory.

If you want to write a script to do this operation (or use ctest),
you have two command-line arguments that overwrite these two options:
- `ReferenceDirectory`
- `OverwriteReference`

Example usage:
`./QAMA -A,MaxIteration=10,ReferenceDirectory="/tmp/my_tmp_refs",OverwriteReference=false ./src/qama/data/tests/data/ExampleFull.arc`


## Load balancing

You can change the loadBalancingMat option (possible values: 'true' / 'false')
and the loadBalancingLoop option (possible values: integer >= 0).

QAMA has two load balancing modes: init (controlled by the loadBalancingMat option)
and compute-loop (controlled by the loadBalancingLoop option).

The first mode is executed once and is only based on material properties.
The second mode is executed every loadBalancingLoop loops (executed if the current
iteration is divisible by loadBalancingLoop) and is based on the number of segments in one loop.


## Input .arc data available

All arc input files examples:

```sh
# (for unbalanced tests, see data/Unbal/README.md)

QS_ARC=${QS_SOURCE_DIR}/data/tests/ExampleFull.arc
QS_ARC=${QS_SOURCE_DIR}/data/tests/ExampleDefault.arc

QS_ARC=${QS_SOURCE_DIR}/data/Strong/Strong_n128.arc
QS_ARC=${QS_SOURCE_DIR}/data/Strong/Strong_n256.arc
QS_ARC=${QS_SOURCE_DIR}/data/Strong/Strong_n512.arc
QS_ARC=${QS_SOURCE_DIR}/data/Strong/Strong_n1024.arc
QS_ARC=${QS_SOURCE_DIR}/data/Strong/Strong_n2048.arc
QS_ARC=${QS_SOURCE_DIR}/data/Strong/Strong_n2520.arc

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