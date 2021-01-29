#!/bin/bash
#PBS -N test
##PBS -l nodes=node02:ppn=64
#PBS -l nodes=1:ppn=8
#PBS -j oe
#PBS -l walltime=5:00:00


export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR
#/opt/software/openmpi-1.10.0/bin/mpirun -np 8 /home/lvp/software/vasp541-me/vasp.5.4.1/bin/vasp_std >log
#mpirun  -np 8 /home/lvp/software/vasp544/opencell/vasp.5.4.4/bin/vasp_std >log
#mpirun  -np 8 /home/lvp/software/vasp544/new-opencell/vasp.5.4.4/bin/vasp_std >log
mpirun  -np 8 /opt/software/vasp544-strain/vasp.5.4.4/bin/vasp_std >log


