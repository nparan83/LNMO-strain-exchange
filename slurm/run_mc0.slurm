#!/bin/bash
#SBATCH -J mc0
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=02:00:00
#SBATCH --output=mc_eps0_%j.out
module purge
module load apps/python/3.11 2>/dev/null || true
t0=$SECONDS
python3 mc.py 0 > mc_eps0.dat
echo "sure: $((SECONDS-t0)) s"
tail -5 mc_eps0.dat
