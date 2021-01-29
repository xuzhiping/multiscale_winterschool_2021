%Frank-Read source

global Bscrew Bedge Beclimb Bline

rn    = [ 1000 1000 1000  7;
		 -1000  -1000 -1000  7;
			0    0     0  0];
links = [1 3 0.5 0.5 0.5 -1 1 0;
		 3 2 0.5 0.5 0.5 -1 1 0];

MU = 1;
NU = 0.305;
maxconnections=8;
lmax = 1000;
lmin = 200;
areamin=lmin*lmin*sin(60/180*pi)*0.5; 
areamax=20*areamin;
a=lmin/sqrt(3)*0.5;
Ec = MU/(4*pi)*log(a/0.1);
totalsteps=200;
dt0=1e7;  
mobility='mobbcc0';
%Drag (Mobility) parameters
Bscrew=1e0;
Bedge=1e0;
Beclimb=1e2;
Bline=1.0e-4*min(Bscrew,Bedge);

integrator='int_trapezoid';
rann = 0.5*a;       
rntol = 0.5*rann;    
doremesh=1;
docollision=1;
doseparation=1;
plotfreq=5;       
plim=10000;         
appliedstress =1e-3.*[2 0 1; 0 2 -1; 1 -1 0];
viewangle=[45 -45 ];
printfreq=1;      
printnode=3;
rmax=100;

% simulation box size (for paradis)
L = 40000;
minSideX = -L/2;
minSideY = -L/2;
minSideZ = -L/2;
maxSideX =  L/2;
maxSideY =  L/2;
maxSideZ =  L/2;

% boundary conditions (for paradis)
xBoundType = 1; % free
yBoundType = 1; % free
zBoundType = 1; % free

% output and communication settings
paradis_dir = '../..';
paradis_input_dir = strcat(paradis_dir,'/Runs');
paradis_output_dir = strcat('Outputs/frank_read_results'); 
