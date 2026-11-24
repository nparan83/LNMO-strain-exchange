#!/bin/bash
#SBATCH -J LNMO_U
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=1-00:00:00
#SBATCH --array=1-8%4
#SBATCH --output=slurm_U_%A_%a.out

export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0

JOB=$(sed -n "${SLURM_ARRAY_TASK_ID}p" U_jobs.txt)
echo "Basliyor: $JOB $(date)"
mpirun pw.x -npool 4 -in ${JOB}.in > ${JOB}.out
RC=$?
echo "Bitti: $JOB rc=$RC $(date)"
SAFE=$(echo $JOB | tr '+-' 'pm')
if [ $RC -eq 0 ]; then rm -rf ./tmp_${SAFE}; fi
