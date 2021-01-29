%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: read_node_data.m
% written by: Seokwoo Lee and Wei Cai
%
% Purpose:  Write ParaDis input files into DDLAB input format
%           (The converted format includes not only nodal information 
%            but also other all parameters.)
%
% How to use: 
%    1. Copy ParaDis file (.data and .cn) to current directory.
%    2. Just change the 'datafile' and 'cnfile' in line 21 and 22.
% 
% Note: 
%    1. Currently, this file can read the fixed format.
%    2. The fixed format is based on the format of output file of ParaDis.
%
%
% Bug_correction
%    1. First version is written - July 17 12:35am, 2007 
%    2. Cosecutively removed nodes (on ParaDis) are corrected - Sep. 22. 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;

% 1. Assign the namve of data and cn file %%%%%%%%%%%%%%%%%%%
datafile = 'rs0001.data';
cnfile = 'restart.cn';


% 2. Generate 'rn matrix' from data file%%%%%%%%%%%%%%%%%%%%%

fid = fopen(datafile,'r');

newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
version_number = fscanf (fid, '%d\n');
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
filesegments_number = fscanf (fid, '%d\n');
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
BoundMin = fscanf (fid, '%e\n');
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
BoundMax = fscanf (fid, '%e\n');
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
Nnodes = fscanf (fid, '%d\n');
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
BoundMinMax = fscanf (fid, '   %e\n');
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
%newline = fgets(fid);

seq = 0;

for i=1:inf
    newline = fscanf (fid, '%e,%e  %e %e %e %e %e \n');
    
    if length(newline) == 15
        seq = seq + 1;
    elseif length(newline) == 8
        seq = seq;
    end
   
    if length(newline) == 15
        rn(seq,:) = [newline(3) newline(4) newline(5) newline(7)];
       
    elseif length(newline) == 0 
       break
    end
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% 2. Generate 'links matrix' from data file%%%%%%%%%%%%%%%%

fid = fopen(datafile,'r');

newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
version_number = fscanf (fid, '%d\n');
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
filesegments_number = fscanf (fid, '%d\n');
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
BoundMin = fscanf (fid, '%e\n');
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
BoundMax = fscanf (fid, '%e\n');
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
Nnodes = fscanf (fid, '%d\n');
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
BoundMinMax = fscanf (fid, '   %e\n');
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
%newline = fgets(fid);

sequ = 0;
seque = 0;


% 2.1 transformation from ParaDis format into DDLAB format

for i=1:inf
    newline = fscanf (fid, '%e,%e  %e %e %e %e %e \n');
        
    sequ = sequ + 1;
   
    if length(newline) == 15
        totallink(sequ,:) = [newline(2)+1 newline(9)+1 ...
            newline(10) newline(11) newline(12) newline(13) newline(14) newline(15)];

    elseif length(newline) == 8
        totallink(sequ,:) = [totallink(sequ-1,1) newline(2)+1 ...
            newline(3) newline(4) newline(5) newline(6) newline(7) newline(8)];
        
    elseif length(newline) == 0 
        break
    end
end

fclose(fid);


% 2.2 Corrections of annihilated nodes in ParaDis
% Sometimes the nodes, which is set initially, are disappeared 
% during dislocation reactions. So, there can be omitted 
% nodes. Thus, it is needed to rearrange the generated links matrix.


% step1. If the tag of first node is not 1, rearrange the totallink matrix 
%        such that it start from 1.
%        Here, totalink matrix is based on ParaDis, so it includes twice
%        number of nodes of DDLAB.

if totallink(1,1)~=1
    decre=totallink(1,1);
    for tt=1:length(totallink)
        totallink(tt,1)=totallink(tt,1)-decre+1;
    end
    
    for jj=1:length(totallink)
        totallink(jj,2)=totallink(jj,2)-decre+1;
    end
end

% step2. Identity the omitted (removed) nodes.
%        Here, 'omit matrix' incldues the tag of the removed nodes.

sequence = 0;
omitseq = 1;
omit = [];

for l=1:length(totallink)
    sequence = sequence + 1;
    if totallink(l,1) == sequence
       sequence = sequence - 1;
    elseif totallink(l,1) ~= sequence
       sequence = sequence + 1;
       if totallink(l,1) == sequence
           sequence = sequence - 1;
       elseif totallink(l,1) ~= sequence
           omit(omitseq,1) = sequence;
           omitseq = omitseq +1;
           if totallink(l,1)-totallink(l-1,1) > 1
               for kk=1:totallink(l,1)-totallink(l-1,1)-2
                   sequence = sequence + 1;
                   if totallink(l,1) ~= sequence;
                       omit(omitseq,1) = sequence;
                       omitseq = omitseq + 1;
                   elseif totallink(l,1) == sequence;
                       sequence = sequence - 1;
                   end
                end
           end
       end
    end
end

% step3. Rearrange the totallink matrix from 1 to end without omiited 
%        nodes.

for rrr=1:length(omit)

            for k=1:length(totallink);
                if totallink(k,1) >= omit(1,1);
                    totallink(k,1) = totallink(k,1) - 1;
                end
            end
                
            for k=1:length(totallink);
                if totallink(k,2) >= omit(1,1);
                    totallink(k,2) = totallink(k,2) - 1;
                end
            end

             
    sequence = 0;
    omitseq = 1;
    omit = [];
    
    for l=1:length(totallink)
        sequence = sequence + 1;
        if totallink(l,1) == sequence
           sequence = sequence - 1;
        elseif totallink(l,1) ~= sequence
           sequence = sequence + 1;
           if totallink(l,1) == sequence
               sequence = sequence - 1;
           elseif totallink(l,1) ~= sequence
               omit(omitseq,1) = sequence;
               omitseq = omitseq +1;
               if totallink(l,1)-totallink(l-1,1) > 1
                   for kk=1:totallink(l,1)-totallink(l-1,1)-2
                       sequence = sequence + 1;
                       if totallink(l,1) ~= sequence;
                           omit(omitseq,1) = sequence;
                           omitseq = omitseq + 1;
                       elseif totallink(l,1) == sequence;
                           sequence = sequence - 1;
                       end
                   end
               end
           end
        end
    end
    
end


% step4. Reduce the matrix
% In ParaDis, each nodal position is appeard twice, but in DDLAB, it
% should be appeared once. Thus, half of nodal positions has to be deleted.

sequen = 1;
for p=1:length(totallink)
    for q=1:length(totallink)
        if [totallink(p,1) totallink(p,2)] == [totallink(q,2) totallink(q,1)];
            Remove(sequen,1) = q;
            sequen = sequen + 1;
        else
        end
    end
end

sequenc = 1;

for h=1:length(Remove)
    if Remove(h,1) > h
        Live(sequenc,1) = h;
        links(sequenc,:) = totallink(h,:);
        sequenc = sequenc + 1;
    end
end


% 3. Generate other needed parameters from cn file%%%%%%%%%
fid = fopen(cnfile,'r');

newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);

maxconnections=8;

lmin = sscanf (fgets(fid), '%*s %*s %e');
lmax = sscanf (fgets(fid), '%*s %*s %e');
a=lmin/sqrt(6);
areamin=lmin*lmin*sin(60/180*pi)*0.5; 
areamax=20*areamin;

dt0=1e-5;

newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);

rmax = sscanf (fgets(fid), '%*s %*s %e');
plotfreq=1;       
plim=5000;

newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
totalsteps = sscanf (fgets(fid), '%*s %*s %d');
MU = sscanf (fgets(fid), '%*s %*s %e');
NU = sscanf (fgets(fid), '%*s %*s %e');
Ec = MU/(4*pi)*log(a/0.1);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);

newline = sscanf (fgets(fid), '%s');
if newline == 'mobilityLaw="FCC_0"'
    mobility = 'mobfcc0';  
elseif newline == 'mobilityLaw="FCC_1"'
   mobility = 'mobfcc1';  
elseif newline == 'mobilityLaw="BCC_0"'
   mobility = 'mobbcc0';  
end

newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);
newline = fgets(fid);

appliedstress = zeros(3,3);

newline = fgets(fid);

appstress = fscanf(fid, '%e');
appliedstress(1,1) = appstress(1);
appliedstress(2,2) = appstress(2);
appliedstress(3,3) = appstress(3);
appliedstress(2,3) = appstress(4);
appliedstress(1,3) = appstress(5);
appliedstress(1,2) = appstress(6);
appliedstress(3,2) = appstress(4);
appliedstress(3,1) = appstress(5);
appliedstress(2,1) = appstress(6);


fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% 4. Another needed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

viewangle=[60 30];
printfreq=1;      
printnode=3;

integrator='int_trapezoid';

rann = 10;       
rntol=1e-1;
%rntol = 2*rann; 

doremesh    =1;
docollision =1;
doseparation=1;

% 5. Variables for ParaDiS
paradis_input_dir = '.';
paradis_output_dir = 'Tests/matlab_results';
minSideX = -17500;
minSideY = -17500;
minSideZ = -17500;
maxSideX = 17500;
maxSideY = 17500;
maxSideZ = 17500;
L = maxSideX - minSideX;
xBoundType = 1;
yBoundType = 1;
zBoundType = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%