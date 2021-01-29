#!/bin/sh
#SBATCH -J GSF_Cu
#SBATCH -p course
#SBATCH -N 1
#SBATCH -o stdout.%j
#SBATCH -e stderr.%j
#SBATCH --no-requeue
#SBATCH --ntasks-per-node=4

module load compiles/intel/2019/u4/config
module load lib/gcc/9.2.0/config
mpiexec.hydra -n 4 /home/train1/WORK/package/lammps-stable_29Oct2020/src/lmp_mpi  < gsf.in > gsf_Cu.log
