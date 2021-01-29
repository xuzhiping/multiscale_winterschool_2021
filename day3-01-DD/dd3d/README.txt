**********************************************************
1. Test cases:

  For example, you can type:
  
    clear all
    run('Inputs\input_frank_read.m');
    dd3d


  More test cases are available in the Inputs/ directory.

**********************************************************
2. Write ParaDiS input files

    clear all
    run('Inputs\input_frank_read.m');
    write_paradis_input

   This will create two files: matlab_input.cn and matlab_input.data
    in the ../../Runs directory, as specified at the bottom of the input file (Inputs\input_frank_read.m).

   This is assuming that the DDLab source code is in the Matlab/dd3d subdirectory of your ParaDiS code.
   Make sure that your ParaDiS directory contains Runs/ and Outputs/ directories.  If not, you can
   create them with the 'mkdir' command.

   You can then run ParaDiS with these input files (outside Matlab):

    cd ~/Codes/ParaDiS
    ./dd3d Runs/matlab_input.cn

   ParaDiS should produce similar results as DDLab.
    
**********************************************************
3. Read ParaDiS output files

   The above ParaDiS run should create a new files in the Outputs/frank_read_results directory,
    such as  Outputs/frank_read_result/restart/restart.cn

   These files can be loaded into DDLab by typing the following in Matlab

    read_node_data('../../Outputs/frank_read_results/restart/restart.data')
    
   [Sorry this one is not working yet, we need to wait for Seokwoo to come back.]
   

**********************************************************
3. You may read about the data structure of the DDLab code from the header of source file

    dd3d.m

