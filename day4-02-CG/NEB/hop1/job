#!/bin/sh
#SBATCH -J lammps-NEB1
#SBATCH -p course
#SBATCH -N 1
#SBATCH -o stdout.%j
#SBATCH -e stderr.%j
#SBATCH --no-requeue
#SBATCH --ntasks-per-node=24

module load compiles/intel/2019/u4/config
module load lib/gcc/9.2.0/config
#mpiexec.hydra -np 24 /home/train1/WORK/package/lammps-stable_29Oct2020/src/lmp_mpi -partition 12x2 -in in.neb.hop1 >log
mpirun -n 24 /home/train1/WORK/package/lammps-stable_29Oct2020/src/lmp_mpi -partition 12x2 -in in.neb.hop1 >log
