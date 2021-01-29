#!/bin/bash
#SBATCH -J 'test'
#SBATCH -p course
#SBATCH -N 1
#SBATCH -o stdout.%J
#SBATCH -e stderr.%J
#SBATCH --no-requeue
#SBATCH --ntasks-per-node=28
module load compiles/intel/2019/u4/config
/apps/soft/abaqus_6.14-4/license/lmgrd -c /apps/soft/abaqus_6.14-4/license/ABAQUS.lic
rm -f abaqus_v6.env
#construct the hosts file for the job

j=''
k=0
scontrol show hostname $SLURM_JOB_NODELIST > ./hosts
line="mp_host_list=["
for i in `cat hosts`
do
  ((k = k + 1))
   if [ $k -ne 1 ]
     then
      line=$line",['"
   else
      line=$line"['"
   fi
   line=$line$i
   line=$line"',28]"
done
line=$line']'
echo $line

sed -n '1,137p' /apps/soft/abaqus_6.14-4/6.14-4/SMA/site/abaqus_v6.env >> abaqus_v6.env
echo $line >> abaqus_v6.env
sed -n '138,$p' /apps/soft/abaqus_6.14-4/6.14-4/SMA/site/abaqus_v6.env >> abaqus_v6.env

/apps/soft/abaqus_6.14-4/Commands/abaqus job=SingleCrack_UEL user=SingleCrack_UEL_PS.f  cpus=28 memory=80000mb interactive
mkdir Result
/apps/soft/abaqus_6.14-4/Commands/abaqus python odb2vtk.py
/apps/soft/abaqus_6.14-4/license/lmdown -c /apps/soft/abaqus_6.14-4/license/ABAQUS.lic -q
