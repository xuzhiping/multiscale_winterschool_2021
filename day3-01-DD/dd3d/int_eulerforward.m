function [rnnew,vn,dt,fn,fseg]=int_eulerforward(rn,dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility)
t=0;


%mobility function
rnvec=[rn(:,1);rn(:,2);rn(:,3)]; flag=rn(:,4);
[vnvec,fn,fseg] = drndt(t,rnvec,flag,MU,NU,a,Ec,links,connectivity,appliedstress,mobility);
%fixed nodes (flag~=0) are not allowed to move
vn=[reshape(vnvec,length(vnvec)/3,3)];
vn=vn.*((rn(:,4)==0)*[1 1 1]);
vmag2=sum(vn.*vn,2);
%limit time step if velocity is too large
vmax=sqrt(max(vmag2));
dt=1/max(1/dt0,vmax/rmax);
%dt=dt0;
rnnew=rn;        
%forward Euler integration
rnnew(:,1:3)=rn(:,1:3)+vn*dt;
    
