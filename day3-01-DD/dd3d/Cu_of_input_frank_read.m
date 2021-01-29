% Frank-Read source

global Bscrew Bedge Beclimb Bline
B = 0.2556e-9;
Bm = 1;%the unit of the length b(amplitude of Burgers vector). 
%When b =10^-9 m,we need to increscent the distance of the immobile node.
%% The definition of the position of initial mobile/mmobile node 
%The size of rn is four columns wide and the number of nodes long.
%The first three columns, contains the x,y,z coordinates of the node, 
%and the fourth column contains a flag. Currently there are only two node flags used in the code.
%A flag equal to zero means that the node is regular node, a flag equal to 7 means that the node is immobile (fixed).
 
%node        x     y    z   flag
% rn    = [ 1000 1000 1000  7;
% 		  -1000  -1000 -1000  7;
% 			0    0     0  0]; 
rn    = [ 3*1*Bm 0*Bm 0*Bm  7;
		 -2*1*Bm  0*Bm 0*Bm  7;
		  0.0*Bm    0.0*Bm     0*Bm  0;
           2*1*Bm 0*Bm 4*1*Bm  7;
		 -3*1*Bm  0*Bm 4*1*Bm  7;
		  0.0*Bm    0.0*Bm     50*Bm  0
           5*1*Bm 0*Bm -4*1*Bm  7;
		 -2*1*Bm  0*Bm -4*1*Bm  7;
		  0.0*Bm    0.0*Bm     -5*1*Bm  0
          2*1*Bm 0*Bm 8*1*Bm  7;
		 -2*1*Bm  0*Bm 8*1*Bm  7;
		  0.0*Bm    0.0*Bm     8*Bm  0
          -4*1*Bm -4*1*Bm -8*1*Bm  7;
		 5*1*Bm  5*1*Bm -8*1*Bm  7;
		  0.0*Bm    0.0*Bm     -8*1*Bm  0];%7 means that the node is immobile

%% The links data structure is eight columns wide and the total number of links long. 
%The first two columns give the node-ids of the starting and ending nodes of the dislocation segment. 
%The 3rd-5th columns of links give the Burgers vector of the dislocation line in Cartesian coordinates,
%and the 6th-8th columns of links gives glide plane of the dislocation segment.

% links = [1 3 0.5 0.5 0.5 -1 1 0;
% 		   3 2 0.5 0.5 0.5 -1 1 0];
%    start_node end_node     b       glide plane
links = [1 3 0.5*Bm 0.5*Bm 0*Bm 0 0 1;
		 3 2 0.5*Bm 0.5*Bm 0*Bm 0 0 1
         4 6 0.5*Bm 0.5*Bm 0*Bm 1 0 1;
		 6 5 0.5*Bm 0.5*Bm 0*Bm 1 0 1
         7 9 0.5*Bm 0.5*Bm 0*Bm 1 0 0;
		 9 8 0.5*Bm 0.5*Bm 0*Bm 1 0 0
         10 12 0.5*Bm 0.5*Bm 0*Bm 1 0 0;
		 12 11 0.5*Bm 0.5*Bm 0*Bm 1 0 0
         13 15 0.5*Bm 0.5*Bm 0.5*Bm 1 0 0;
		 15 14 0.5*Bm 0.5*Bm 0.5*Bm 1 0 0];    
%%  materials parameters    
MU = 48e9;% shear model(Pa)
NU = 0.3;%  Poisson's ratio
Y_e = 2*MU*(1+NU);% Young¡¯s modulus(Pa)
maxconnections=8; % maximum number of segments a node may have
lmax = 2*Bm; %for remesh,maximum length of a dislocation segment(b) 
lmin = 1*Bm;  %for remesh,minimum length of a dislocation segment (b) 
areamin=lmin*lmin*sin(60/180*pi)*0.5; % minimum area criterion,for remesh
areamax=20*areamin; 
a=0.1*Bm; %lmin/sqrt(3)*0.5;  % dislocation core radius used for non-singular force calculation
Ec = MU/(4*pi)*log(0.1/0.1); % dislocation core energy per unit length
totalsteps=200; % number of cycles that are run for completion of dd3d command
dt0= 2*10^-10; % time steep (s)
mobility='mobbcc0'; %mobility law
%Drag (Mobility) parameters
Bscrew=1;%3.3333*10^-7Pa.s; %Drag coefficient
Bedge=1;%3.3333*10^-7Pa.s; %Drag coefficient
Beclimb=1e0;
Bline=1.0e0*min(Bscrew,Bedge);
rann = 0.5*a;  % annihilation distance used to calculate the collision of dislocation lines     
rntol = 0.5*rann;   
%% Input for the topology
doremesh=1; % do remesh
docollision=1; %do collision
doseparation=1;%do separation  ;deal with junction?
integrator='int_trapezoid';% Selection of numerical methods
%% output 
plotfreq=5; % number of cycles between plots of geometry      
plim=20*Bm; % limits of plotting space (box of the  crystal)    
viewangle=[45 -45 ];
printfreq=1;  % number of cycles between plots of geometry    
printnode=3;  % nodeid of the monitored node
rmax=100*Bm; 
%% appliedstress
%appliedstress0 =1e-3.*[2 0 1; 0 2 -1; 1 -1 0]; 
appliedstress0 =10^-1*[0 0 0; 0 0 0; 0 0 0]*MU; %initial stress field (Pa)
stress_rate = 1*10^7*MU;% stress rate
%% simulation box size (for paradis)
L = 4*Bm;%10^-6m
minSideX = -L/2;
minSideY = -L/2;
minSideZ = -L/2;
maxSideX =  L/2;
maxSideY =  L/2;
maxSideZ =  L/2;
%% boundary conditions (for paradis)
xBoundType = 1; % free
yBoundType = 1; % free
zBoundType = 1; % free
%% output and communication settings
paradis_dir = '../..';
paradis_input_dir = strcat(paradis_dir,'/Runs');
paradis_output_dir = strcat('Outputs/frank_read_results'); 
