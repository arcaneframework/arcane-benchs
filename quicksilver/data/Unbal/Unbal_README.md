```
UnbalX_Y_nZ
X = Resolution
Y = num of cellules / sub-domain (csd)
Z = num of procs

1000 particles/cell

1csd =
  1 cell/sd (1x1x1)
  1000 particles/sd

16csd =
  16 cells/sd (4x2x2)
  16000 particles/sd

(
  Example: Unbal32_16_n32768:
    X = 32
    Y = 16 cells/sd (4x2x2)
    Z = 32768

    32/4 = 8

    32768 = 8*64*64

    num sd x: 8
    num sd y: 64
    num sd z: 64

    num cells x: 8*4 = 32
    num cells y: 64*2 = 128
    num cells z: 64*2 = 128
)


H = Hard (split++)
S = Soft (absorb++)

----

    / H / S / H / S /|
   / H / S / H / S / |
  / H / S / H / S /  |
 / H / S / H / S /   |
| H | S | H | S |   /
| H | S | H | S |  /
| H | S | H | S | /
| H | S | H | S |/

z y
|/
--x
(4x4x4 cells)

----
Resolution 4:

      | H | S | H | S |
       -> x
1csd: |000|001|002|003|

Unbal4_1:
  length(x) = 4*4 = 16
  length(y) = length(x)*8 = 128
  length(z) = length(x)*8 = 128
Unbal4_1_n4.arc:     nSub-domains: 4x1x1     nCells: 4x1x1    nParticles: 4000
Unbal4_1_n16.arc:    nSub-domains: 4x2x2     nCells: 4x2x2    nParticles: 16000
Unbal4_1_n32.arc:    nSub-domains: 4x4x2     nCells: 4x4x2    nParticles: 32000
Unbal4_1_n128.arc:   nSub-domains: 4x8x4     nCells: 4x8x4    nParticles: 128000
Unbal4_1_n512.arc:   nSub-domains: 4x16x8    nCells: 4x16x8   nParticles: 512000
Unbal4_1_n2048.arc:  nSub-domains: 4x32x16   nCells: 4x32x16  nParticles: 2048000
Unbal4_1_n16384.arc: nSub-domains: 4x64x64   nCells: 4x64x64  nParticles: 16384000
Unbal4_1_n32768.arc: nSub-domains: 4x128x64  nCells: 4x1x1    nParticles: 32768000

----

----
Resolution 16:

       | S | H | H | S | S | S | S | H | H | H | H | S | S | S | S | S |
        -> x
1csd:  |000|001|002|003|004|005|006|007|008|009|010|011|012|013|014|015|
16csd: |      000      |      001      |      002      |      003      |

Unbal16:
  length(x) = 16*4 = 64
  length(y) = length(x)*8 = 512
  length(z) = length(x)*8 = 512
Unbal16_1_n16.arc:     nSub-domains: 16x1x1    nCells: 16x1x1    nParticles: 16000
Unbal16_1_n32.arc:     nSub-domains: 16x2x1    nCells: 16x2x1    nParticles: 32000
Unbal16_1_n128.arc:    nSub-domains: 16x4x2    nCells: 16x4x2    nParticles: 128000
Unbal16_1_n512.arc:    nSub-domains: 16x8x4    nCells: 16x8x4    nParticles: 512000
Unbal16_1_n2048.arc:   nSub-domains: 16x16x8   nCells: 16x16x8   nParticles: 2048000
Unbal16_1_n16384.arc:  nSub-domains: 16x32x32  nCells: 16x32x32  nParticles: 16384000
Unbal16_1_n32768.arc:  nSub-domains: 16x64x32  nCells: 16x64x32  nParticles: 32768000

Unbal16_16_n4.arc:     nSub-domains: 4x1x1     nCells: 16x2x2      nParticles: 64000
Unbal16_16_n16.arc:    nSub-domains: 4x2x2     nCells: 16x4x4      nParticles: 256000
Unbal16_16_n32.arc:    nSub-domains: 4x4x2     nCells: 16x8x4      nParticles: 512000
Unbal16_16_n128.arc:   nSub-domains: 4x8x2     nCells: 16x16x8     nParticles: 2048000
Unbal16_16_n512.arc:   nSub-domains: 4x16x8    nCells: 16x32x16    nParticles: 8192000
Unbal16_16_n2048.arc:  nSub-domains: 4x32x16   nCells: 16x64x32    nParticles: 32768000
Unbal16_16_n16384.arc: nSub-domains: 4x64x64   nCells: 16x128x128  nParticles: 262144000
Unbal16_16_n32768.arc: nSub-domains: 4x128x64  nCells: 16x256x128  nParticles: 524288000

----

----
Resolution 32:

       | S | H | H | S | H | S | S | S | S | S | S | S | S | S | H | H | H | H | H | H | S | H | S | S | S | S | S | S | H | S | H | S |
        -> x
1csd:  |000|001|002|003|004|005|006|007|008|009|010|011|012|013|014|015|016|017|018|019|020|021|022|023|024|025|026|027|028|029|030|031|
16csd: |      000      |      001      |      002      |      003      |      004      |      005      |      006      |      007      |

Unbal32:
  length(x) = 32*4 = 128
  length(y) = length(x)*8 = 1024
  length(z) = length(x)*8 = 1024
Unbal32_1_n32.arc:     nSub-domains: 32x1x1    nCells: 32x1x1    nParticles: 32000
Unbal32_1_n128.arc:    nSub-domains: 32x2x2    nCells: 32x2x2    nParticles: 128000
Unbal32_1_n512.arc:    nSub-domains: 32x4x4    nCells: 32x4x4    nParticles: 512000
Unbal32_1_n2048.arc:   nSub-domains: 32x8x8    nCells: 32x8x8    nParticles: 2048000
Unbal32_1_n16384.arc:  nSub-domains: 32x32x16  nCells: 32x32x16  nParticles: 16384000
Unbal32_1_n32768.arc:  nSub-domains: 32x32x32  nCells: 32x32x32  nParticles: 32768000

Unbal32_16_n8.arc:     nSub-domains: 8x1x1    nCells: 32x2x2      nParticles: 128000
Unbal32_16_n32.arc:    nSub-domains: 8x2x2    nCells: 32x4x4      nParticles: 512000
Unbal32_16_n128.arc:   nSub-domains: 8x4x4    nCells: 32x8x8      nParticles: 2048000
Unbal32_16_n512.arc:   nSub-domains: 8x8x8    nCells: 32x16x16    nParticles: 8192000
Unbal32_16_n2048.arc:  nSub-domains: 8x16x16  nCells: 32x32x32    nParticles: 32768000
Unbal32_16_n16384.arc: nSub-domains: 8x64x32  nCells: 32x128x64   nParticles: 262144000
Unbal32_16_n32768.arc: nSub-domains: 8x64x64  nCells: 32x128x128  nParticles: 524288000

----

mpirun -n 4 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal4_1/Unbal4_1_n4.arc
mpirun -n 16 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal4_1/Unbal4_1_n16.arc
mpirun -n 32 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal4_1/Unbal4_1_n32.arc
mpirun -n 128 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal4_1/Unbal4_1_n128.arc
mpirun -n 512 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal4_1/Unbal4_1_n512.arc
mpirun -n 2048 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal4_1/Unbal4_1_n2048.arc
mpirun -n 16384 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal4_1/Unbal4_1_n16384.arc
mpirun -n 32768 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal4_1/Unbal4_1_n32768.arc

mpirun -n 16 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_1/Unbal16_1_n16.arc
mpirun -n 32 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_1/Unbal16_1_n32.arc
mpirun -n 128 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_1/Unbal16_1_n128.arc
mpirun -n 512 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_1/Unbal16_1_n512.arc
mpirun -n 2048 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_1/Unbal16_1_n2048.arc
mpirun -n 16384 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_1/Unbal16_1_n16384.arc
mpirun -n 32768 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_1/Unbal16_1_n32768.arc

mpirun -n 4 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_16/Unbal16_16_n4.arc
mpirun -n 16 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_16/Unbal16_16_n16.arc
mpirun -n 32 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_16/Unbal16_16_n32.arc
mpirun -n 128 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_16/Unbal16_16_n128.arc
mpirun -n 512 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_16/Unbal16_16_n512.arc
mpirun -n 2048 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_16/Unbal16_16_n2048.arc
mpirun -n 16384 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_16/Unbal16_16_n16384.arc
mpirun -n 32768 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal16_16/Unbal16_16_n32768.arc

mpirun -n 32 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal32_1/Unbal32_1_n32.arc
mpirun -n 128 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal32_1/Unbal32_1_n128.arc
mpirun -n 512 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal32_1/Unbal32_1_n512.arc
mpirun -n 2048 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal32_1/Unbal32_1_n2048.arc
mpirun -n 16384 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal32_1/Unbal32_1_n16384.arc
mpirun -n 32768 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal32_1/Unbal32_1_n32768.arc

mpirun -n 4 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal32_16/Unbal32_16_n8.arc
mpirun -n 32 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal32_16/Unbal32_16_n32.arc
mpirun -n 128 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal32_16/Unbal32_16_n128.arc
mpirun -n 512 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal32_16/Unbal32_16_n512.arc
mpirun -n 2048 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal32_16/Unbal32_16_n2048.arc
mpirun -n 16384 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal32_16/Unbal32_16_n16384.arc
mpirun -n 32768 ${QS_EXE} ${QS_SOURCE_DIR}/data/Unbal/Unbal32_16/Unbal32_16_n32768.arc
```
