#!/bin/bash
#BATCH -J vasp                     # job name 
#SBATCH -p course                  # queue  --cnall,  course
#SBATCH -N 1                       # number of nodes requested
#SBATCH -n 28                       # 
#SBATCH -o stdout.%j
#SBATCH -e stderr.%j
#SBATCH --no-requeue
module load compiles/intel/2019/u4/config
mpiexec.hydra -n 28 /apps/soft/vasp/vasp.5.4.4/e5_2680v4/opa/vasp.5.4.4/bin/vasp_std > log
