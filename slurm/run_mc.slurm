#!/bin/bash
#SBATCH -J lnmo_mc
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=04:00:00
#SBATCH --array=0-10
#SBATCH --output=slurm_mc_%a_%A.out
E=(-5 -4 -3 -2 -1 0 1 2 3 4 5)
X=${E[$SLURM_ARRAY_TASK_ID]}
t0=$SECONDS
python3 mc2.py $X > mc_eps${X}.dat
echo "eps=$X sure: $((SECONDS-t0)) s"
