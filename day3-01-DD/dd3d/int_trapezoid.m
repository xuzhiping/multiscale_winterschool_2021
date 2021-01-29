function [rn,vn,dt,fn,fseg]=int_trapezoid(rn,dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility)
%implicit numerical integrator using the Trapezoid method
%dt: suggested timestep (usually from previous iteration)
%dt0: maximum allowed timestep

%dummy variable
t=0;
rnold=rn;

%scramble rn into a single column vector
rnvec=[rn(:,1);rn(:,2);rn(:,3)]; flag=rn(:,4);

%Backward Euler
rnvec0=rnvec;

[vnvec0,fn,fseg]=drndt(t,rnvec0,flag,MU,NU,a,Ec,links,connectivity,appliedstress,mobility);
%dt=1/max(1/dt0,max(vnvec0)/rmax);
%dt=dt0;

dt1=dt;
maxiter=1;
convergent=0;
while(~convergent)
    rnvec1=rnvec0+vnvec0*dt;

    for iter=1:maxiter,
        [vnvec,fn,fseg]=drndt(t,rnvec1,flag,MU,NU,a,Ec,links,connectivity,appliedstress,mobility);
        %err=rnvec1-rnvec0-vnvec.*dt;          %backward Euler
        err=rnvec1-rnvec0-(vnvec+vnvec0)/2*dt; %trapzoid
        errmag=max(abs(err));
        disp(sprintf('iter=%d err=%e',iter,errmag));
        if(errmag<rntol)
            convergent=1;
            break;
        else
            rnvec1=rnvec1-err;
        end
    end
    
    if(convergent)
        break;
    else
        dt=dt/2;
    end
end

%unscramble rn and vn vectors
rn=[reshape(rnvec1,length(rnvec1)/3,3),flag];
vn=reshape(vnvec,length(vnvec)/3,3); % trapezoidal rule modification
    
%automatically adjust time step
if((dt==dt1)&(iter==1))
    maxchange=1.2;
    exponent=20;
    factor=maxchange*(1/(1+(maxchange^exponent-1)*(errmag/rntol)))^(1/exponent);
    dt=min(dt1*factor,dt0);
end