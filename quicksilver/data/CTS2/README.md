CTS2

1 rank/core
4096 cells/rank
40960 particles/rank

-----------------
-----------------

QS_EXE= #QS arcane executable
QS_SOURCE_DIR= #QS arcane source path

ORI_QS_EXE= #QS original executable
ORI_QS_SOURCE_DIR= #QS original source path

unset OMP_NUM_THREADS

mpirun -n 128    ${QS_EXE} ${QS_SOURCE_DIR}/data/CTS2/CTS2_n128.arc
mpirun -n 128    ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CTS2_Benchmark/CTS2.inp  -X 128 -Y 64 -Z 64     -x 128 -y 64 -z 64     -I 8 -J 4 -K 4     -n 5242880

mpirun -n 256    ${QS_EXE} ${QS_SOURCE_DIR}/data/CTS2/CTS2_n256.arc
mpirun -n 256    ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CTS2_Benchmark/CTS2.inp  -X 128 -Y 128 -Z 64    -x 128 -y 128 -z 64    -I 8 -J 8 -K 4     -n 10485760

mpirun -n 512    ${QS_EXE} ${QS_SOURCE_DIR}/data/CTS2/CTS2_n512.arc
mpirun -n 512    ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CTS2_Benchmark/CTS2.inp  -X 128 -Y 128 -Z 128   -x 128 -y 128 -z 128   -I 8 -J 8 -K 8     -n 20971520

mpirun -n 1024   ${QS_EXE} ${QS_SOURCE_DIR}/data/CTS2/CTS2_n1024.arc
mpirun -n 1024   ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CTS2_Benchmark/CTS2.inp  -X 256 -Y 128 -Z 128   -x 256 -y 128 -z 128   -I 16 -J 8 -K 8    -n 41943040

mpirun -n 2048   ${QS_EXE} ${QS_SOURCE_DIR}/data/CTS2/CTS2_n2048.arc
mpirun -n 2048   ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CTS2_Benchmark/CTS2.inp  -X 256 -Y 256 -Z 128   -x 256 -y 256 -z 128   -I 16 -J 16 -K 8   -n 83886080

mpirun -n 4096   ${QS_EXE} ${QS_SOURCE_DIR}/data/CTS2/CTS2_n4096.arc
mpirun -n 4096   ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CTS2_Benchmark/CTS2.inp  -X 256 -Y 256 -Z 256   -x 256 -y 256 -z 256   -I 16 -J 16 -K 16  -n 167772160

mpirun -n 8192   ${QS_EXE} ${QS_SOURCE_DIR}/data/CTS2/CTS2_n8192.arc
mpirun -n 8192   ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CTS2_Benchmark/CTS2.inp  -X 512 -Y 256 -Z 256   -x 512 -y 256 -z 256   -I 32 -J 16 -K 16  -n 335544320

mpirun -n 131072 ${QS_EXE} ${QS_SOURCE_DIR}/data/CTS2/CTS2_n131072.arc
mpirun -n 131072 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CTS2_Benchmark/CTS2.inp  -X 2048 -Y 512 -Z 512  -x 2048 -y 512 -z 512  -I 64 -J 64 -K 32  -n 5368709120
