#!/bin/bash
#SBATCH -J lnmo_cfg0
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=08:00:00
#SBATCH --output=slurm_cfg0_%j.out

export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0

mkdir -p out
t0=$SECONDS
mpirun pw.x -npool 4 -in eps0_cfg0.in > eps0_cfg0.out
echo "cfg0 sure: $((SECONDS-t0)) s"
grep '!    total energy'      eps0_cfg0.out | tail -1
grep 'total magnetization'    eps0_cfg0.out | tail -1
grep 'absolute magnetization' eps0_cfg0.out | tail -1
