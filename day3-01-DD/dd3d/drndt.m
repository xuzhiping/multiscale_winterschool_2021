function [vnvec,fn,fseg] = drndt(t,rnvec,flag,MU,NU,a,Ec,links,connectivity,appliedstress,mobility)

%unscramble rn
rn=[reshape(rnvec,length(rnvec)/3,3),flag];

%rn(:,1:3)

%nodal driving force
fseg=segforcevec(MU,NU,a,Ec,rn,links,appliedstress,0);
    
%mobility function
[vn,fn]=feval(mobility,fseg,rn,links,connectivity,[],[]);

%fixed nodes (flag~=0) are not allowed to move
vn=vn.*((rn(:,4)==0)*[1 1 1]);    

%make up the last
vnvec=reshape(vn,length(vn(:,1))*3,1);

%vn
%disp(sprintf('t=%20.10e vn(1,:)=(%20.10e,%20.10e,%20.10e)',t,vn(1,1),vn(1,2),vn(1,3)));
%pause(0.2)
%if(t>0)
%   pause
%end