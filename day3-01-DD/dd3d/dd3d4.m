%% Dislocation Dynamics simulation in 3-dimension
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
%% compile the c source code for seg-seg force evaluation and makes a dynamic linked library
mex -O SegSegForces.c
 
%default value if run by itself (e.g. not through "rundd3d")
%% cleanup the empty node and link entries at the end of the initial data structures
[rn,links]=cleanupnodes(rn,links);%
%% genererate the connectivity list from the list of links
[connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections); 
consistencycheck(rn,links,connectivity,linksinconnect);
%% plot dislocation structure
figure(1); 
subplot(1,2,1)
plotnodes(rn,links,plim); view(viewangle); xlim([-plim plim]); ylim([-plim plim]); zlim([-plim plim]);
title('The dislocation configuration');
subplot(1,2,2)
plot ([0,0],[0,0],'Color','r','MarkerSize',25);
hold on;
xlim([0 0.5]); ylim([0 0.5]); 
xlabel('\epsilon_{xz}/1');ylabel('\tau_{xz}/MU')
title('The stress-stain curve');
drawnow
frame=getframe(gcf);
im{1}=frame2im(frame);
if ~exist('./../Result','dir')
    mkdir ./../Result
end
if ~exist('./../Result/example-4','dir')
    mkdir ./../Result/example-4
end

data=zeros(totalsteps,1);
%% step time 
if(~exist('dt'))
    dt=dt0;
end
dt=min(dt,dt0);


mdold=10;
counter=0;
%% variables for recording the stress and strain
density =zeros(totalsteps,1);
appliedstress=appliedstress0;
stress11 = zeros (totalsteps,1);
stress12 = zeros (totalsteps,1);
stress13 = zeros (totalsteps,1);
stress22 = zeros (totalsteps,1);
stress23 = zeros (totalsteps,1);
stress33 = zeros (totalsteps,1);
strain11 = zeros (totalsteps,1);
strain12 = zeros (totalsteps,1);
strain13 = zeros (totalsteps,1);
strain22 = zeros (totalsteps,1);
strain23 = zeros (totalsteps,1);
strain33 = zeros (totalsteps,1);

ep11 = zeros (totalsteps,1);
ep_11 = 0;
ep12 = zeros (totalsteps,1);
ep_12 = 0;
ep13 = zeros (totalsteps,1);
ep_13 = 0;
ep22 = zeros (totalsteps,1);
ep_22 = 0;
ep23 = zeros (totalsteps,1);
ep_23 = 0;
ep33 = zeros (totalsteps,1);
ep_33 = 0;

strain_11 = 0;
strain_12 = 0;
strain_13 = 0;
strain_22 = 0;
strain_23 = 0;
strain_33 = 0;

dstress11 = zeros (totalsteps,1);
dstress12 = zeros (totalsteps,1);
dstress13 = zeros (totalsteps,1);
dstress22 = zeros (totalsteps,1);
dstress23 = zeros (totalsteps,1);
dstress33 = zeros (totalsteps,1);

sum11 = zeros (totalsteps,1);
sum12 = zeros (totalsteps,1);
sum13 = zeros (totalsteps,1);
sum22 = zeros (totalsteps,1);
sum23 = zeros (totalsteps,1);
sum33 = zeros (totalsteps,1);
%e = zeros (totalsteps,1);
%delta_t = zeros(totalsteps,1);
%% cycle of the all step time
for curstep=1:totalsteps  
    %% integrating equation of motion
    [rnnew,vn,dt,fn,fseg]=feval(integrator,rn,dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility);%计算速度
    % plastic strain and plastic spin calculations  
    %% increment of the stress in dt
    dstress = dt*(stress_rate)*[0 0 1; 0 0 0; 1 0 0];
    appliedstress =appliedstress+dstress;
    %delta_t(curstep+1,1) = dt;    
    %% record of stress and strain 
    stress11(curstep+1,1) = appliedstress(1,1);
    stress12(curstep+1,1) = appliedstress(1,2);
    stress13(curstep+1,1) = appliedstress(1,3);
    stress22(curstep+1,1) = appliedstress(2,2);
    stress23(curstep+1,1) = appliedstress(2,3);
    stress33(curstep+1,1) = appliedstress(3,3);
    
    dstress11(curstep+1,1) = dstress(1,1);
    dstress12(curstep+1,1) = dstress(1,2);
    dstress13(curstep+1,1) = dstress(1,3);
    dstress22(curstep+1,1) = dstress(2,2);
    dstress23(curstep+1,1) = dstress(2,3);
    dstress33(curstep+1,1) = dstress(3,3);
    
    strain11(curstep+1,1) = (1/Y_e)*dstress11(curstep+1,1)-(NU/Y_e)*dstress22(curstep+1,1)-(NU/Y_e)*dstress33(curstep+1,1);
    strain22(curstep+1,1) = -(NU/Y_e)*dstress11(curstep+1,1)+(1/Y_e)*dstress22(curstep+1,1)-(NU/Y_e)*dstress33(curstep+1,1);
    strain33(curstep+1,1) = -(NU/Y_e)*dstress11(curstep+1,1)-(NU/Y_e)*dstress22(curstep+1,1)+(1/Y_e)*dstress33(curstep+1,1);
    strain12(curstep+1,1) = (1+NU)/Y_e*dstress12(curstep+1,1);
    strain13(curstep+1,1) = (1+NU)/Y_e*dstress13(curstep+1,1)/1;
    strain23(curstep+1,1) = (1+NU)/Y_e*dstress23(curstep+1,1);
    
    strain_11 = strain_11 + strain11(curstep+1,1);
    sum11(curstep+1,1) = strain_11;
    strain_12 = strain_12 + strain12(curstep+1,1);
    sum12(curstep+1,1) = strain_12;
    strain_13 = strain_13 + strain13(curstep+1,1);
    sum13(curstep+1,1) = strain_13;
    strain_22 = strain_22 + strain13(curstep+1,1);
    sum13(curstep+1,1) = strain_22;
    strain_23 = strain_23 + strain23(curstep+1,1);
    sum23(curstep+1,1) = strain_23;
    strain_33 = strain_33 + strain33(curstep+1,1);
    sum33(curstep+1,1) = strain_33;
    %% plastic strain and plastic spin calculations
    [ep_inc,wp_inc]=calcplasticstrainincrement(rnnew,rn,links,(2*plim)^3);%计算应变
    
    ep11(curstep+1,1) = ep_inc(1,1);
    ep_11 = ep_11+ep_inc(1,1);
    ep12(curstep+1,1) = ep_inc(1,2);
    ep_12 = ep_12+ep_inc(1,2);    
    ep13(curstep+1,1) = ep_inc(1,3);
    ep_13 = ep_13+ep_inc(1,3);
    ep22(curstep+1,1) = ep_inc(2,2);
    ep_22 = ep_13+ep_inc(2,2);
    ep23(curstep+1,1) = ep_inc(2,3);
    ep_23 = ep_13+ep_inc(2,3);
    ep33(curstep+1,1) = ep_inc(3,3);
    ep_33 = ep_13+ep_inc(3,3);
    %% plot stress_strain curve
    for ii = curstep-1:curstep
    if(mod(curstep,plotfreq)==0)
      subplot(1,2,2)
      plot ([0,sum13(curstep+1,1)-ep13(curstep+1,1)],[0,appliedstress(1,3)/(MU)],'Color','r','MarkerSize',25);
      hold on;
      xlim([0 0.5]); ylim([0 0.5]); 
      xlabel('\epsilon_{xz}/1');ylabel('\tau_{xz}/MU')
      title('The stress-stain curve');
      drawnow

      pause(0.04);
    end
    end
    %% print step,dt and the velocity of node
    if(mod(curstep,printfreq)==0)
        disp(sprintf('step%3d dt=%e v%d=(%e,%e,%e)',...
            curstep,dt,printnode,vn(printnode,1),vn(printnode,2),vn(printnode,3)));%显示每个变量的值
    end
    %% plot the evolution of the dislocation loop
    if(mod(curstep,plotfreq)==0)
        subplot(1,2,1)
        plotnodes(rn,links,plim);  xlim([-plim plim]); ylim([-plim plim]); zlim([-plim plim]);
        view(viewangle);
        title('The dislocation configuration');
        drawnow
        frame=getframe(gcf);
        im{round(curstep/plotfreq+1)}=frame2im(frame);       
        pause(0.01);
    end
    %% for remseh (update the node position)
    rnnew=[rnnew(:,1:3) vn rnnew(:,4)];%多加入了速度
    linksnew=links;
    connectivitynew=connectivity;
    linksinconnectnew=linksinconnect;
    fsegnew=fseg;
    if(doseparation)
        %spliting of nodes with 4 or more connections
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=separation(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,mobility,MU,NU,a,Ec,2*rann,appliedstress);
    %将连接多的点分成两个
    end
    if(docollision)
        %collision detection and handling
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=collision(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,appliedstress);
    end
    
    if(doremesh)
        %remesh
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=remesh(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,appliedstress,mobility);      
    end
   
    rn=[rnnew(:,1:3) rnnew(:,7)];%多加了速度
    vn=rnnew(:,4:6);
    links=linksnew;
    connectivity=connectivitynew;
    linksinconnect=linksinconnectnew;
    %% Calculate the dislocation density
    LL = 0;
    for i = 4:length(links(:,1))
        %zeros(length(links(:,1)),1);
        ll = rn(links(i,2),1:3)-rn(links(i,1),1:3);
        LL = sqrt(dot(ll,ll))+LL;       
    end
    V = plim*plim*plim;
    density(curstep,1) = LL/V;
    
    fseg=fsegnew;
    %% store run time information
    %time step
    data(curstep,1)=dt;
    save restart
end
save restart
for idx = 1:size(im,2)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,['./../Result/example-4/conf.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
    else
        imwrite(A,map,['./../Result/example-4/conf.gif'],'gif','WriteMode','append','DelayTime',0.5);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
