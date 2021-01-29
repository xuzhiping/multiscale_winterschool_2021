#!/bin/sh
#SBATCH -J cnt_tear
#SBATCH -p course
#SBATCH -N 1
#SBATCH -o stdout.%j
#SBATCH -e stderr.%j
#SBATCH --no-requeue
#SBATCH --ntasks-per-node=28

module load compiles/intel/2019/u4/config
#exe="/apps/soft/lammps/lammps-7Aug19/e5_2680v4/opa/lammps-7Aug19/src/lmp_mpi"
#exe="/apps/soft/lammps/lammps-3Mar20/e5_2680v4/opa/lammps-7Aug19/src/lmp_mpi"
exe="/home/train1/WORK/package/lammps-stable_29Oct2020/src/lmp_mpi"

mpiexec.hydra -n 24 ${exe} < cnt_tear.in >& log
