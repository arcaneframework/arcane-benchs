Coral2_P1

4096 cells/node
40 particles/cell

-----------------
Coral2_P2

1331 cells/node
40 particles/cell

-----------------
-----------------

Coral2_P1

QS_EXE= #QS arcane executable
QS_SOURCE_DIR= #QS arcane source path

ORI_QS_EXE= #QS original executable
ORI_QS_SOURCE_DIR= #QS original source path

export OMP_NUM_THREADS=2
mpirun -n 128 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n128.arc
mpirun -n 128 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem1/Coral2_P1.inp  -X 16 -Y 16 -Z 16  -x 16 -y 16 -z 16  -I 8 -J 4 -K 4  -n 163840

mpirun -n 256 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n256.arc
mpirun -n 256 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem1/Coral2_P1.inp  -X 32 -Y 16 -Z 16  -x 32 -y 16 -z 16  -I 8 -J 8 -K 4  -n 327680

mpirun -n 512 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n512.arc
mpirun -n 512 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem1/Coral2_P1.inp  -X 32 -Y 32 -Z 16  -x 32 -y 32 -z 16  -I 8 -J 8 -K 8  -n 655360

mpirun -n 1024 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n1024.arc
mpirun -n 1024 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem1/Coral2_P1.inp  -X 32 -Y 32 -Z 32  -x 32 -y 32 -z 32  -I 16 -J 8 -K 8  -n 1310720

mpirun -n 2048 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n2048.arc
mpirun -n 2048 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem1/Coral2_P1.inp  -X 64 -Y 32 -Z 32  -x 64 -y 32 -z 32  -I 16 -J 16 -K 8  -n 2621440

mpirun -n 4096 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n4096.arc
mpirun -n 4096 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem1/Coral2_P1.inp  -X 64 -Y 64 -Z 32  -x 64 -y 64 -z 32  -I 16 -J 16 -K 16  -n 5242880

mpirun -n 8192 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n8192.arc
mpirun -n 8192 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem1/Coral2_P1.inp  -X 64 -Y 64 -Z 64  -x 64 -y 64 -z 64  -I 32 -J 16 -K 16  -n 10485760

mpirun -n 131072 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P1_n131072.arc
mpirun -n 131072 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem1/Coral2_P1.inp  -X 512 -Y 256 -Z 256  -x 512 -y 256 -z 256  -I 64 -J 64 -K 32  -n 167772160


Coral2_P2

export OMP_NUM_THREADS=2
mpirun -n 128 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n128.arc
mpirun -n 128 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem2/Coral2_P2.inp  -X 11 -Y 11 -Z 11  -x 11 -y 11 -z 11  -I 8 -J 4 -K 4  -n 53240

mpirun -n 256 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n256.arc
mpirun -n 256 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem2/Coral2_P2.inp  -X 22 -Y 11 -Z 11  -x 22 -y 11 -z 11  -I 8 -J 8 -K 4  -n 106480

mpirun -n 512 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n512.arc
mpirun -n 512 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem2/Coral2_P2.inp  -X 22 -Y 22 -Z 11  -x 22 -y 22 -z 11  -I 8 -J 8 -K 8  -n 212960

mpirun -n 1024 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n1024.arc
mpirun -n 1024 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem2/Coral2_P2.inp  -X 22 -Y 22 -Z 22  -x 22 -y 22 -z 22  -I 16 -J 8 -K 8  -n 425920

mpirun -n 2048 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n2048.arc
mpirun -n 2048 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem2/Coral2_P2.inp  -X 44 -Y 22 -Z 22  -x 44 -y 22 -z 22  -I 16 -J 16 -K 8  -n 851840

mpirun -n 4096 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n4096.arc
mpirun -n 4096 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem2/Coral2_P2.inp  -X 44 -Y 44 -Z 22  -x 44 -y 44 -z 22  -I 16 -J 16 -K 16  -n 1703680

mpirun -n 8192 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n8192.arc
mpirun -n 8192 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem2/Coral2_P2.inp  -X 44 -Y 44 -Z 44  -x 44 -y 44 -z 44  -I 32 -J 16 -K 16  -n 3407360

mpirun -n 131072 ${QS_EXE} -A,T=2 ${QS_SOURCE_DIR}/data/Coral2/Coral2_P2_n131072.arc
mpirun -n 131072 ${ORI_QS_EXE} -i ${ORI_QS_SOURCE_DIR}/Examples/CORAL2_Benchmark/Problem2/Coral2_P2.inp  -X 176 -Y 88 -Z 88  -x 176 -y 88 -z 88  -I 64 -J 64 -K 32  -n 6814720
