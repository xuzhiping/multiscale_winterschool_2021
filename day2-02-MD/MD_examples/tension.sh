#!/bin/sh
#SBATCH -J MD_tension
#SBATCH -p course
#SBATCH -N 1
#SBATCH -o stdout.%j
#SBATCH -e stderr.%j
#SBATCH --no-requeue
#SBATCH --ntasks-per-node=4

module load compiles/intel/2019/u4/config
module load lib/gcc/9.2.0/config

scontrol show hostname $SLURM_JOB_NODELIST > ./hosts
rm -rf hostfile
touch hostfile
for i in `cat hosts`
do
 echo $i:4>>hostfile
done

mpirun -machinefile hostfile -np 4 /home/train1/WORK/package/lammps-stable_29Oct2020/src/lmp_mpi  < tension.in > tension_Cu.log
