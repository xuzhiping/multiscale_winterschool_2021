#!/bin/bash
#SBATCH -J 'damask'
#SBATCH -p course
#SBATCH -N 1
#SBATCH -o stdout.%J
#SBATCH -e stderr.%J
#SBATCH --no-requeue
#SBATCH --ntasks-per-node=6

#---------------------------------------------------------------
## Step 1: Set the execution environment

source /apps/soft/DAMASK-v2.0.1/DAMASK_env.sh
export PETSC_DIR=/apps/src/petsc-3.6.4/
export PETSC_ARCH=arch-linux2-c-opt
export PATH=/apps/src/petsc-3.6.4/arch-linux2-c-opt/bin:$PATH
export PATH=/apps/lib/anaconda/anaconda3/e5/bin:$PATH
export DAMASK_NUM_THREADS=6

#---------------------------------------------------------------
## Step 2: Run the spectral solver for calculation

/apps/soft/DAMASK-v2.0.1/code/DAMASK_spectral.exe --load tensionX.load --geom Myself.geom

#---------------------------------------------------------------
## Step 3: Postprocessing the data of various fields

module load lib/anaconda/anaconda2/config

# Output the data
postResults Myself_tensionX.spectralOut --increments --range 0 400 100 --separation x,y,z --split --cr texture,f,p,orientation,grainrotation

cd postProc/

# Processing data at increment 0
addStrainTensors Myself_tensionX_inc000.txt --left --logarithmic;
addCauchy Myself_tensionX_inc000.txt;
vtk_rectilinearGrid Myself_tensionX_inc000.txt;
vtk_addRectilinearGridData -s 1_Cauchy,2_Cauchy,3_Cauchy,4_Cauchy,5_Cauchy,6_Cauchy,1_orientation,2_orientation,3_orientation,1_grainrotation,2_grainrotation,3_grainrotation Myself_tensionX_inc000.txt --vtk Myself_tensionX_inc000_pos\(cell\).vtr;

# Processing data at increment 100
addStrainTensors Myself_tensionX_inc100.txt --left --logarithmic;
addCauchy Myself_tensionX_inc050.txt;
vtk_rectilinearGrid Myself_tensionX_inc100.txt;
vtk_addRectilinearGridData -s 1_Cauchy,2_Cauchy,3_Cauchy,4_Cauchy,5_Cauchy,6_Cauchy,1_orientation,2_orientation,3_orientation,1_grainrotation,2_grainrotation,3_grainrotation Myself_tensionX_inc100.txt --vtk Myself_tensionX_inc100_pos\(cell\).vtr;

# Processing data at increment 200
addStrainTensors Myself_tensionX_inc200.txt --left --logarithmic;
addCauchy Myself_tensionX_inc200.txt;
vtk_rectilinearGrid Myself_tensionX_inc200.txt;
vtk_addRectilinearGridData -s 1_Cauchy,2_Cauchy,3_Cauchy,4_Cauchy,5_Cauchy,6_Cauchy,1_orientation,2_orientation,3_orientation,1_grainrotation,2_grainrotation,3_grainrotation Myself_tensionX_inc200.txt --vtk Myself_tensionX_inc200_pos\(cell\).vtr;


# Processing data at increment 300
addStrainTensors Myself_tensionX_inc300.txt --left --logarithmic;
addCauchy Myself_tensionX_inc300.txt;
vtk_rectilinearGrid Myself_tensionX_inc300.txt;
vtk_addRectilinearGridData -s 1_Cauchy,2_Cauchy,3_Cauchy,4_Cauchy,5_Cauchy,6_Cauchy,1_orientation,2_orientation,3_orientation,1_grainrotation,2_grainrotation,3_grainrotation Myself_tensionX_inc300.txt --vtk Myself_tensionX_inc300_pos\(cell\).vtr;

# Processing data at increment 400
addStrainTensors Myself_tensionX_inc400.txt --left --logarithmic;
addCauchy Myself_tensionX_inc400.txt;
vtk_rectilinearGrid Myself_tensionX_inc400.txt;
vtk_addRectilinearGridData -s 1_Cauchy,2_Cauchy,3_Cauchy,4_Cauchy,5_Cauchy,6_Cauchy,1_orientation,2_orientation,3_orientation,1_grainrotation,2_grainrotation,3_grainrotation Myself_tensionX_inc400.txt --vtk Myself_tensionX_inc400_pos\(cell\).vtr;

cd ..

mv postProc FieldData


#---------------------------------------------------------------
## Step 4: Postprocessing the data vs time

postResults Myself_tensionX.spectralOut --cr f,p --co edge_density --co dipole_density --co twin_fraction

cd postProc

addStrainTensors Myself_tensionX.txt --left --logarithmic

addCauchy Myself_tensionX.txt

addMises -e 'ln(V)' -s Cauchy Myself_tensionX.txt

cd ..

mv postProc HistoryData
