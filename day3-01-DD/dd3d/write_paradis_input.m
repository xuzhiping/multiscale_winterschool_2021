%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: write_paradis_input.m
%
% Purpose:  Write simulation setting into ParaDiS format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% write text input (restart) file for paradis
fname = strcat(paradis_input_dir,'/matlab_input.cn');
fp = fopen(fname,'w');

fprintf(fp,'###ParaDiS restart file created by DDLab (use -*-shell-script-*- format)\n');

fprintf(fp,'\n');
fprintf(fp,'#Directory to write output files\n');
fprintf(fp,'dirname = %s\n',paradis_output_dir);

fprintf(fp,'\n');
fprintf(fp,'#Input files for PBC image stresses\n');
fprintf(fp,'Rijmfile = "Inputs/Rijm.cube.out"\n');
fprintf(fp,'RijmPBCfile = "Inputs/RijmPBC.cube.out"\n');

fprintf(fp,'\n');
fprintf(fp,'#Total simulation steps\n');
fprintf(fp,'maxstep = %d\n',totalsteps);

fprintf(fp,'\n');
fprintf(fp,'#The total number of CPUs should be numXdoms * numYdoms * numZdoms\n');
fprintf(fp,'#numXdoms = 1\n');
fprintf(fp,'#numYdoms = 1\n');
fprintf(fp,'#numZdoms = 1\n');

fprintf(fp,'\n');
fprintf(fp,'#Cells for dislocation grouping (cannot be less than 3)\n');
fprintf(fp,'numXcells = 3\n');
fprintf(fp,'numYcells = 3\n');
fprintf(fp,'numZcells = 3\n');

fprintf(fp,'\n');
fprintf(fp,'#Fast multipole method specs.\n');
fprintf(fp,'fmEnabled = 0  #disable fast multipole\n');
fprintf(fp,'fmMPOrder = 2\n');
fprintf(fp,'fmTaylorOrder = 4\n');
fprintf(fp,'fmCorrectionTbl = "Inputs/fm-ctab.m2.t4.dat"\n');

fprintf(fp,'\n');
fprintf(fp,'timestepIntegrator = "backward-euler"\n');

fprintf(fp,'\n');
fprintf(fp,'#Simulation box size\n');
fprintf(fp,'xBoundMin = %f\n',minSideX);
fprintf(fp,'xBoundMax = %f\n',maxSideX);
fprintf(fp,'yBoundMin = %f\n',minSideY);
fprintf(fp,'yBoundMax = %f\n',maxSideY);
fprintf(fp,'zBoundMin = %f\n',minSideZ);
fprintf(fp,'zBoundMax = %f\n',maxSideZ);

fprintf(fp,'\n');
fprintf(fp,'#Boundary conditions\n');
fprintf(fp,'xBoundType = %d\n',xBoundType);
fprintf(fp,'yBoundType = %d\n',yBoundType);
fprintf(fp,'zBoundType = %d\n',zBoundType);

fprintf(fp,'\n');
fprintf(fp,'#Fundamental length unit\n');
fprintf(fp,'burgMag = 1.0\n');

fprintf(fp,'\n');
fprintf(fp,'#Elastic constants\n');
fprintf(fp,'shearModulus = %f\n',MU);
fprintf(fp,'pois = %f\n',NU);

fprintf(fp,'\n');
fprintf(fp,'#Mobility law function\n');
if (strcmp(mobility,'mobfcc0'))
    fprintf(fp,'mobilityLaw = "%s"\n','FCC_0');
elseif (strcmp(mobility,'mobbcc0'))
    fprintf(fp,'mobilityLaw = "%s"\n','BCC_0');
else
    disp(sprintf('unknown mobility law (%s)! use FCC_0 in ParaDiS',mobility));
    fprintf(fp,'mobilityLaw = "%s"\n','FCC_0');
end

fprintf(fp,'\n');
fprintf(fp,'MobScrew = %.6e\n',1/Bscrew);
fprintf(fp,'MobEdge  = %.6e\n',1/Bedge);
fprintf(fp,'MobClimb = %.6e\n',1/Beclimb);

fprintf(fp,'\n');
fprintf(fp,'#Discretization\n');
fprintf(fp,'maxSeg = %f\n',lmax);
fprintf(fp,'minSeg = %f\n',lmin);

fprintf(fp,'\n');
fprintf(fp,'#Maximum nodal displacement at each time step\n');
fprintf(fp,'rmax = %g\n',rmax);

fprintf(fp,'\n');
fprintf(fp,'#Error tolerance in determining time step\n');
fprintf(fp,'rTol = %g\n',rntol);

fprintf(fp,'\n');
fprintf(fp,'#Maximum time step\n');
fprintf(fp,'maxDT = %g\n',dt0);

fprintf(fp,'\n');
fprintf(fp,'#Core cut-off radius\n');
fprintf(fp,'rc = %f\n\n',a);

fprintf(fp,'\n');
fprintf(fp,'#Core energy\n');
fprintf(fp,'Ecore = %e\n',Ec);

fprintf(fp,'\n');
fprintf(fp,'#Turn on elastic interaction\n');
fprintf(fp,'elasticinteraction = 1\n');

fprintf(fp,'\n');
fprintf(fp,'#Applied stress in Pa (xx,yy,zz,yz,zx,xy)\n');
fprintf(fp,'appliedStress = [ %g %g %g %g %g %g ]\n', ...
            appliedstress(1), appliedstress(2), appliedstress(3), ...
            appliedstress(4), appliedstress(5), appliedstress(6) );

fprintf(fp,'\n');
fprintf(fp,'#Save files\n');
fprintf(fp,'savecn = 1\n');
fprintf(fp,'savecnfreq = 100\n');
            
fprintf(fp,'\n');
fprintf(fp,'preserveOldTags = 1\n');
%fprintf(fp,'command_file = "command_pipe"\n');


fclose(fp);

disp(sprintf('ParaDiS input file %s created.\n',fname));

% create node data file
fname = strcat(paradis_input_dir,'/matlab_input.data');
write_node_data(fname,rn,links,L);
disp(sprintf('ParaDiS node data file %s created.\n',fname));


