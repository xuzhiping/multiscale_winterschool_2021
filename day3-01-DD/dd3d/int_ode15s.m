function [rn,vn,dt,fn,fseg]=int_ode15s(rn,dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility)

%dummy variable
t=0;
rnold=rn;

%scramble rn into a single column vector
rnvec=[rn(:,1);rn(:,2);rn(:,3)]; flag=rn(:,4);

if(0) %Forward Euler
    %time derivative
    [vnvec,fn,fseg]=drndt(t,rnvec,flag,MU,NU,a,Ec,links,connectivity,appliedstress,mobility);
    %limit time step if velocity is too large
    vmax=max(max(abs(vnvec)));
    dt=1/max(1/dt0,vmax/rntol);
    %forward Euler integration
    rnvec=rnvec+vnvec*dt;
end

%this method is inefficient because it uses finite difference to find Jacobian (second derivative)
options=odeset('RelTol',rntol);
[T,solution]=ode15s(@drndt,[0 dt0],rnvec,options,flag,MU,NU,a,links,connectivity,appliedstress,mobility);
%[T,solution]=ode23t(@drndt,[0 dt0/100],rnvec,options,flag,MU,NU,a,links,connectivity,nsub,appliedstress,mobility);
%[T,solution]=ode23(@drndt,[0 dt0/10],rnvec,[],flag,MU,NU,a,links,connectivity,nsub,appliedstress,mobility);

T,size(solution)

dt=T(end);
rnvec=solution(end,:)';

%unscramble rn and vn vectors
rn=[reshape(rnvec,length(rnvec)/3,3),flag];
vn=(rn(:,1:3)-rnold(:,1:3))/dt;
    
