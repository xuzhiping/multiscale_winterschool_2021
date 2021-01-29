%Multi Junction Formation

global Bscrew Bedge Beclimb Bline

disp(' ');
disp(' In this example, the multi-arm node is supposed to split by itself');
disp(' Unfortunately, this is not working yet.  The split nodes automatically recombines.');
disp(' We will debug this test case in the future...');

theta = 35.2644*pi/180;
theta= 90*pi/180;
L = 1000;
vz = 1/sqrt(3) .* [ 1  1  1];
v1 = 1/sqrt(6) .* [-2  1  1];
v2 = 1/sqrt(6) .* [ 1 -2  1];
v3 = 1/sqrt(6) .* [ 1  1 -2];
l1 = L .* (cos(theta) .* vz + sin(theta) .* v1);
l2 = L .* (cos(theta) .* vz + sin(theta) .* v2);
l3 = L .* (cos(theta) .* vz + sin(theta) .* v3);
links = [1 7 -0.5 0.5 0.5 0 1 -1;
         7 2 -0.5 0.5 0.5 0 1 -1;
         3 7 0.5 -0.5 0.5 1 0 -1;
         7 4 0.5 -0.5 0.5 1 0 -1;
         5 7 0.5 0.5 -0.5 1 -1 0;
         7 6 0.5 0.5 -0.5 1 -1 0];
rn    = [-0.5 .*l1 7;
          0.5 .*l1 7;
         -0.5 .*l2 7;
          0.5 .*l2 7;
         -0.5 .*l3 7;
          0.5 .*l3 7;
          0  0  0  0];
      
maxconnections=8;
lmax = 600;
lmin = 200;
a=lmin/sqrt(6);
MU = 1.3e11;
NU = 0.309; % for BCC Mo
Ec = MU/(4*pi)*log(a/0.1);
totalsteps=200;
areamin=lmin*lmin*sin(60/180*pi)*0.5; % minimum discretization area
areamax=20*areamin; % maximum discretization area
dt0=200000;           %maximum time step
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
rann = 1e-3;       %annihilation distance (capture radius)
%=1e-1;       %tolerance for integrating equation of motion
rntol = 0.5*rann;      % on Tom's suggestion
doremesh =1;
docollision=1;
doseparation=1;
dt=1e-1;


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