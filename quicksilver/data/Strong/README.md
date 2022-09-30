```
Strong (Scaling)

128 ranks/node
2 threads/rank
262144 cells
128000000 particles


-----------------
-----------------

QS_EXE= #QS arcane executable
QS_SOURCE_DIR= #QS arcane source path

ORI_QS_EXE= #QS original executable
ORI_QS_SOURCE_DIR= #QS original source path

export OMP_NUM_THREADS=2

mpirun -n 128    ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Strong/Strong_n128.arc
mpirun -n 128    ${ORI_QS_EXE} -i ${QS_SOURCE_DIR}/data/Strong/input_for_qs/Strong.inp -I 8 -J 4 -K 4

mpirun -n 256    ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Strong/Strong_n256.arc
mpirun -n 256    ${ORI_QS_EXE} -i ${QS_SOURCE_DIR}/data/Strong/input_for_qs/Strong.inp -I 8 -J 8 -K 4

mpirun -n 512    ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Strong/Strong_n512.arc
mpirun -n 512    ${ORI_QS_EXE} -i ${QS_SOURCE_DIR}/data/Strong/input_for_qs/Strong.inp -I 8 -J 8 -K 8

mpirun -n 1024   ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Strong/Strong_n1024.arc
mpirun -n 1024   ${ORI_QS_EXE} -i ${QS_SOURCE_DIR}/data/Strong/input_for_qs/Strong.inp -I 16 -J 8 -K 8

mpirun -n 2048   ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Strong/Strong_n2048.arc
mpirun -n 2048   ${ORI_QS_EXE} -i ${QS_SOURCE_DIR}/data/Strong/input_for_qs/Strong.inp -I 16 -J 16 -K 8

mpirun -n 2520   ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Strong/Strong_n2520.arc
mpirun -n 2520   ${ORI_QS_EXE} -i ${QS_SOURCE_DIR}/data/Strong/input_for_qs/Strong.inp -I 16 -J 16 -K 10

```