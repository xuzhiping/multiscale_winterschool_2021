#!/bin/bash
#SBATCH -J matlab
#SBATCH -p course
#SBATCH -N 1
#SBATCH -o stdout.%J
#SBATCH -e stderr.%J
#SBATCH --no-requeue
#SBATCH --ntasks-per-node=1
/apps/soft/MATLAB/R2018a/bin/matlab -nodesktop -nosplash -r "run p3.m;quit"

