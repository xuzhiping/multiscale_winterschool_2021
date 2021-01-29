#!/bin/bash
#BATCH -J vasp                     # job name 
#SBATCH -p course                  # queue  --cnall,  course
#SBATCH -N 1                       # number of nodes requested
#SBATCH -n 12                       # 
#SBATCH -o stdout.%j
#SBATCH -e stderr.%j
#SBATCH --no-requeue
module load compiles/intel/2019/u4/config
mpiexec.hydra -n 12 /home/train1/bin/vasp.5.4.4/bin/vasp_std > log
