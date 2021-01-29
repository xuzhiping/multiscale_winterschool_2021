%Binary Junction Formation

global Bscrew Bedge Beclimb Bline

rn=[ 4000          0       4000      7;
     2001.9847     7.0235  2001.9849 0;
    -4000.0000     0.0000 -4000.0000 7;
        0.0000  4000.0000  4000.0000 7;
      -69.4714   -69.4714   -69.4740 0;
        0.0000 -4000.0000 -4000.0000 7;
    -2001.5798    -6.4075 -2001.5799 0;
        7.0233  2001.9848  2001.9850 0;
       -6.4067 -2001.5798 -2001.5800 0;
        75.5022   75.5022    75.5050 0];
  b1=[-1 1 -1]./2;
  b2=[1 -1 -1]./2;
  b3=[0 0 1];
  n1=[1 1 0]./sqrt(2);
  n2=[0 1 1]./sqrt(2);
  n3=[1 0 1]./sqrt(2);
  
links=[1 2   b1 n1;
       2 10  b1 n1;
       5 7   b1 n1;
       7 3   b1 n1;
       4 8   b2 n2;
       8 10  b2 n2;
       5 9   b2 n2;
       9 6   b2 n2;
       5 10  b3 n1];
       
     
maxconnections=8;
lmax = 3000;
lmin = 1000;
a=lmin/sqrt(6);
MU = 1.3e11;
NU = 0.309; % for BCC Mo
Ec = MU/(4*pi)*log(a/0.1);
totalsteps=200;
areamin=lmin*lmin*sin(60/180*pi)*0.5; % minimum discretization area
areamax=20*areamin; % maximum discretization area
dt0=1e-5;           %maximum time step
rmax=10.0;         %maximum allowed displacement per timestep
plotfreq=5;       %plot nodes every how many steps
plim=5000;          %plot x,y,z limit (nodes outside not plotted)
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
%mobility = 'mobbcc1b';
%Drag (Mobility) parameters
Bscrew=1e0;
Bedge=1e0;
Beclimb=1e2;
Bline=1.0e-4*min(Bscrew,Bedge);

integrator='int_trapezoid';
rann = 10;       %annihilation distance (capture radius)
%rntol=1e-1;       %tolerance for integrating equation of motion
rntol = 2*rann;      % on Tom's suggestion
rmax=3;
doremesh =1;
docollision=1;
doseparation=1;
dt=1e-10;

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
