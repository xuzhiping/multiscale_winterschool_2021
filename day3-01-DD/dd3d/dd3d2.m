

%Dislocation Dynamics simulation in 3-dimension
%
%the best way to run it is by "rundd3d"
%
%Features:
%infinite boundary condition (no pbc)
%linear mobility law (mobfcc0,mobfcc1)
%no remesh, no collision detection, no reaction
%N^2 interaction (no neighbor list, no fast multipole)

%Data structure:
%NMAX:    maximum number of nodes (including disabled ones)
%LINKMAX: maximum number of links (including disabled ones)
%rn: (NMAX,4) array of nodal positions (last column is flag: -1 means disabled)
%vn: (NMAX,3) array of nodal velocities
%fn: (NMAX,3) array of nodal forces
%links: (LINKMAX,8) array of links (id1,id2,bx,by,bz,nx,ny,nz)

% compile the c source code for seg-seg force evaluation and makes a dynamic linked library
mex -O SegSegForces.c

%default value if run by itself (e.g. not through "rundd3d")
% cleanup the empty node and link entries at the end of the initial data structures
[rn,links]=cleanupnodes(rn,links);

% genererate the connectivity list from the list of links
[connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
consistencycheck(rn,links,connectivity,linksinconnect);

%plot dislocation structure
figure(1); 
plotnodes(rn,links,plim); view(viewangle); xlim([-plim plim]); ylim([-plim plim]); zlim([-plim plim]);
drawnow
frame=getframe(gcf);
im{1}=frame2im(frame);
if ~exist('./../Result','dir')
    mkdir ./../Result
end
if ~exist('./../Result/example-2','dir')
    mkdir ./../Result/example-2
end
data=zeros(totalsteps,1);
if(~exist('dt'))
    dt=dt0;
end
dt=min(dt,dt0);
mdold=10;
counter=0;
for curstep=1:totalsteps
    %integrating equation of motion
    [rnnew,vn,dt,fn,fseg]=feval(integrator,rn,dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility);
    % plastic strain and plastic spin calculations
    [ep_inc,wp_inc]=calcplasticstrainincrement(rnnew,rn,links,(2*plim)^3);
    
    
    if(mod(curstep,printfreq)==0)
        disp(sprintf('step%3d dt=%e v%d=(%e,%e,%e)',...
            curstep,dt,printnode,vn(printnode,1),vn(printnode,2),vn(printnode,3)));
    end
    if(mod(curstep,plotfreq)==0)
        plotnodes(rn,links,plim);  xlim([-plim plim]); ylim([-plim plim]); zlim([-plim plim]);
        view(viewangle);
        drawnow
        frame=getframe(gcf);
        im{round(curstep/plotfreq+1)}=frame2im(frame);
        pause(0.01);
    end
    
    rnnew=[rnnew(:,1:3) vn rnnew(:,4)];
    linksnew=links;
    connectivitynew=connectivity;
    linksinconnectnew=linksinconnect;
    fsegnew=fseg;
    if(doseparation)
        %spliting of nodes with 4 or more connections
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=separation(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,mobility,MU,NU,a,Ec,2*rann,appliedstress);
    
    end
    if(docollision)
        %collision detection and handling
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=collision(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,appliedstress);
    end
    
    if(doremesh)
        %remesh
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=remesh(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,appliedstress,mobility);      
    end
   
    rn=[rnnew(:,1:3) rnnew(:,7)];
    vn=rnnew(:,4:6);
    links=linksnew;
    connectivity=connectivitynew;
    linksinconnect=linksinconnectnew;
    fseg=fsegnew;
    %store run time information
    %time step
    data(curstep,1)=dt;
    save restart
end
save restart
for idx = 1:size(im,2)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,['./../Result/example-2/conf.gif'],'gif','LoopCount',Inf,'DelayTime',0.01);
    else
        imwrite(A,map,['./../Result/example-2/conf.gif'],'gif','WriteMode','append','DelayTime',0.01);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%