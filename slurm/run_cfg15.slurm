#!/bin/bash
#SBATCH -J lnmo_cfg
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=03:00:00
#SBATCH --array=1-5
#SBATCH --output=slurm_cfg%a_%A.out
export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
C=eps0_cfg${SLURM_ARRAY_TASK_ID}
mkdir -p out
t0=$SECONDS
mpirun pw.x -npool 4 -in $C.in > $C.out
echo "$C sure: $((SECONDS-t0)) s"
grep -E "convergence has been achieved" $C.out | tail -1
grep -E "^!    total energy" $C.out | tail -1
grep -E "total magnetization" $C.out | tail -1
