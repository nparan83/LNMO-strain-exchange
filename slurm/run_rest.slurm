#!/bin/bash
#SBATCH -J lnmo_rest
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=03:00:00
#SBATCH --array=0-31%8
#SBATCH --output=slurm_rest_%a_%A.out
export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
E=(-04 -03 -02 -01 +01 +02 +03 +04)
M=(FM A C G)
i=$SLURM_ARRAY_TASK_ID
C=e${E[$((i/4))]}_${M[$((i%4))]}
mkdir -p out
t0=$SECONDS
mpirun pw.x -npool 4 -in $C.in > $C.out
echo "$C sure: $((SECONDS-t0)) s"
grep -E "convergence has been achieved" $C.out | tail -1
grep -E "^!    total energy" $C.out | tail -1
grep -E "total magnetization" $C.out | tail -1
