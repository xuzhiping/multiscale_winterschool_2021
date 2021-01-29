%Binary Junction Formation

global Bscrew Bedge Beclimb Bline

rn=[ 61 -60 300 7 ;
     21 -20 100 0;
    -19  20 -100 0 ;
    -59  60 -300 7 ;
     60  60 -330 7 ;
     20  20 -110 0;
    -20 -20  110 0;
    -60 -60  330 7] ;
links=[1 2 0.5 -0.5 -0.5 1 1 0;
       2 3 0.5 -0.5 -0.5 1 1 0;
       3 4 0.5 -0.5 -0.5 1 1 0;
       5 6 0.5 0.5 -0.5 1 -1 0;
       6 7 0.5 0.5 -0.5 1 -1 0;
       7 8 0.5 0.5 -0.5 1 -1 0];
   
 rn=[6.100000000000000e+001   -6.000000000000000e+001    3.000000000000000e+002  7;
     5.0e-001    5.0e-001   -2.613534012614e+000  0; 
    -6.000000000000000e+001   -6.000000000000000e+001    3.300000000000000e+002  7;
    -5.900000000000000e+001    6.000000000000000e+001   -3.000000000000000e+002  7;
     6.000000000000000e+001    6.000000000000000e+001   -3.300000000000000e+002  7];
links=[1.0e+000    2.0e+000    5.0e-001;
          2.0e+000    4.0e+000    5.0e-001;
    2.0e+000    3.0e+000    5.0e-001;
    5.0e+000    2.0e+000    5.0e-001];
    
links(:,4:6)=[    -5.0e-001   -5.0e-001    1.0e+000;
   -5.0e-001   -5.0e-001    1.0e+000;
    5.0e-001   -5.0e-001    1.0e+000;
    5.0e-001   -5.0e-001    1.0e+000];
 links(:,7:8)=[    1.0e+000                         0;
    1.0e+000                         0;
   -1.0e+000                         0;
   -1.0e+000                         0];
maxconnections=8;
lmax = 500;
lmin = 100;
a=lmin/sqrt(6);
MU = 1.3e11;
NU = 0.309; % for BCC Mo
Ec = MU/(4*pi)*log(a/0.1);
totalsteps=200;
areamin=lmin*lmin*sin(60/180*pi)*0.5; % minimum discretization area
areamax=20*areamin; % maximum discretization area
dt0=1e-5;           %maximum time step
rmax=10.0;         %maximum allowed displacement per timestep
plotfreq=1;       %plot nodes every how many steps
plim=500;          %plot x,y,z limit (nodes outside not plotted)
appliedstress = zeros(3,3);
%sigma = [ 1 0 0 ; 0 1 0 ; 0 0 -2];
%sigma = [ 0 0 0 ; 0 0 0 ; 0 0  1];
%sigma = [ 0 1 0 ; 1 0 1 ; 0 1 0];
sigma = [ 1 -1 0 ; -1 1 0 ; 0 0 0];
viewangle=[-45 45];
printfreq=1;      %print out information every how many steps
printnode=3;
%mobility='mobfcc0'; %mobility law function
mobility='mobbcc0';
%Drag (Mobility) parameters
Bscrew=1e0;
Bedge=1e0;
Beclimb=1e2;
Bline=1.0e-4*min(Bscrew,Bedge);

integrator='int_trapezoid';
%mobility = 'mobbcc1b';
rann = 2;       %annihilation distance (capture radius)
%rntol=1e-1;       %tolerance for integrating equation of motion
rntol = 1;      % on Tom's suggestion
doremesh =1;
docollision=1;
doseparation=1;
dt=1e-10;

% simulation box size (for paradis)
L = 4000;
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
