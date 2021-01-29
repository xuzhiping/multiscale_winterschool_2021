%Multi Junction Formation

global Bscrew Bedge Beclimb Bline

rn=[-100 -100    0 7;
       0 -100 -100 7;
    -100    0 -100 7;
     100  100    0 7;
       0  100  100 7;
     100    0  100 7;
     -50  -50  -50 0;
      50   50   50 0];
  b1=[1 1 -1]./sqrt(3);
  b2=[-1 1 1]./sqrt(3);
  b3=[1 -1 1]./sqrt(3);
  b4=[1 1 1]./sqrt(3);
  n1=[1 -1 0]./sqrt(2);
  n2=[0 -1 1]./sqrt(2);
  n3=[1 0 -1]./sqrt(2);
  
links=[1 7  b1 n1;
       2 7  b2 n2;
       3 7  b3 n3;
       4 8 -b1 n1;
       5 8 -b2 n2;
       6 8 -b3 n3;
       7 8  b4 n1];
     
maxconnections=8;
lmax = 50;
lmin = 10;
a=lmin/sqrt(6);
MU = 1.3e11;
NU = 0.309; % for BCC Mo
Ec = MU/(4*pi)*log(a/0.1);
totalsteps=10000;   % changed!cyn
areamin=lmin*lmin*sin(60/180*pi)*0.5; % minimum discretization area
areamax=20*areamin; % maximum discretization area
dt0=20000;           %maximum time step
rmax=10.0;         %maximum allowed displacement per timestep
plotfreq=1000;       %plot nodes every how many steps
plim=150;          %plot x,y,z limit (nodes outside not plotted)
%sigma = [ 1 0 0 ; 0 1 0 ; 0 0 -2];
sigma = [ 0 0 0 ; 0 0 0 ; 0 0  1];
%sigma = [ 0 1 0 ; 1 0 1 ; 0 1 0];
%sigma = [ 1 -1 0 ; -1 1 0 ; 0 0 0];
appliedstress=zeros(3,3);
%appliedstress=2.3e9.*sigma;   % changed!cyn
appliedstress=5e8.*sigma;   % changed!cyn
appliedstressrate=zeros(3,3);   % changed!cyn
%appliedstressrate=1e12.*sigma;   % changed!cyn
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
rann = 1.0e-2;       %annihilation distance (capture radius)
%rntol=1e-1;       %tolerance for integrating equation of motion
rntol = rann;      % on Tom's suggestion
doremesh =1;
docollision=1;
doseparation=1;
dt=1e-10;

% simulation box size (for paradis)
L = 2000;
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
