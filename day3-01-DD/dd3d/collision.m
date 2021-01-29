function [rn,links,connectivity,linksinconnect,fseg]=collision(rn,links,connectivity,linksinconnect,fseg,mindist,MU,NU,a,Ec,mobility,appliedstress);
% this subroutine goes through the existing links and checks for collisions
% it first checks through unconnected links
% it then checks for hinges that are coming within the mininum distance
mindist2=mindist*mindist;
lrn2=length(rn(1,:));
lrn3=lrn2-1;
% eps is the error factor of the calculation
eps=1e-12;
% check for two links colliding
i=1;
while i<size(links,1)
    j=i+1;
    while j<=size(links,1)
        n1s1=links(i,1);
        n2s1=links(i,2);
        n1s2=links(j,1);
        n2s2=links(j,2);
        if (n1s1~=n1s2)&(n1s1~=n2s2)&(n2s1~=n1s2)&(n2s1~=n2s2)
            [dist2,ddist2dt,L1,L2]=mindistcalc(rn(n1s1,1:lrn3),rn(n2s1,1:lrn3),rn(n1s2,1:lrn3),rn(n2s2,1:lrn3));
            collision_condition_is_met=((dist2<mindist2)&(ddist2dt<-eps))|(dist2<eps);
           
            % there are two conditions here the first condition handles non planar collisions
            % the second conditions looks for coplanar collisions
            if collision_condition_is_met
                % links are unconnected and colliding
                % identify the first node to be merged
                vec=rn(n1s1,1:3)-rn(n2s1,1:3);
                close_to_n1s1=((L1*L1*(vec*vec'))<mindist2);
                close_to_n2s1=(((1-L1)*(1-L1)*(vec*vec'))<mindist2);
                % if collision point is close to one of the existing nodes use that node
                if close_to_n1s1
                    mergenode1=n1s1;
                elseif close_to_n2s1
                    mergenode1=n2s1;
                else
                    spnode=n1s1;
                    splitconnection=linksinconnect(i,1);
                    posvel=rn(n1s1,1:lrn3).*(1-L1)+rn(n2s1,1:lrn3).*L1;
                    [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel);
                    mergenode1=length(rn(:,1));
                    linknew=length(links(:,1));
                    links(linknew,6:8)=links(i,6:8);
                    fseg=[fseg;zeros(1,6)];
                end

                % identify the second node to be merged
                vec=rn(n1s2,1:3)-rn(n2s2,1:3);
                close_to_n1s2=((L2*L2*(vec*vec'))<mindist2);
                close_to_n2s2=(((1-L2)*(1-L2)*(vec*vec'))<mindist2);
                % if collision point is close to one of the existing nodes use that node
                if close_to_n1s2 
                    mergenode2=n1s2;
                elseif close_to_n2s2
                    mergenode2=n2s2;
                else
                    spnode=n1s2;
                    splitconnection=linksinconnect(j,1);
                    posvel=rn(n1s2,1:lrn3).*(1-L2)+rn(n2s2,1:lrn3).*L2;
                    [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel);
                    mergenode2=length(rn(:,1));
                    linknew=length(links(:,1));
                    links(linknew,6:8)=links(j,6:8);
                    fseg=[fseg;zeros(1,6)];
                end
                % merge the two colliding nodes
                disp(sprintf('node %d and node %d are colliding by two line collision',mergenode1,mergenode2))
                collisionpoint=findcollisionpoint(mergenode1,mergenode2,rn,connectivity,links);
                rn(mergenode1,1:lrn2)=[collisionpoint 0 0 0 max(rn(mergenode1,lrn2),rn(mergenode2,lrn2)) ];
                [rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,mergenode1,mergenode2,MU,NU,a,Ec);
                if mergednodeid>0
                    for k=1:connectivity(mergednodeid,1)
                        linkid=connectivity(mergednodeid,2*k);
                        fseg(linkid,:)=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),links,appliedstress,linkid);
                        othernode=links(linkid,3-connectivity(mergednodeid,2*k+1));
                        clist=[connectivity(othernode,1) linspace(1,connectivity(othernode,1),connectivity(othernode,1))];
                        [rn(othernode,4:6),fntmp]=feval(mobility,fseg,rn,links,connectivity,othernode,clist);
                    end
                    numbcon=connectivity(mergednodeid,1);
                    conlist=[numbcon linspace(1,numbcon,numbcon)];
                    [rn(mergednodeid,4:6),fntmp]=feval(mobility,fseg,rn,links,connectivity,mergednodeid,conlist);
                end
            end
        end
        j=j+1;
    end
    i=i+1;
end
% check for a hinge condition
i=1;
while i<=length(rn(:,1))
    j=1;
    while j<=connectivity(i,1)
        nodenoti=links(connectivity(i,2*j),3-connectivity(i,2*j+1));
        k=1;
        while k<=connectivity(i,1)
            linkid=connectivity(i,2*k);
            % if node is on the link do not check for collision
            if j~=k
                % identify the nodes on the link
                n1s1=links(linkid,1);
                n2s1=links(linkid,2);
                [dist2,ddist2dt,L1,L2]=mindistcalc(rn(n1s1,1:lrn3),rn(n2s1,1:lrn3),rn(nodenoti,1:lrn3),rn(nodenoti,1:lrn3));
                collision_condition_is_met=(dist2<mindist2)&(ddist2dt<-eps);
                if collision_condition_is_met
                    % identify the first node to be merged
                    mergenode1=nodenoti;
                    % identify the second node to be merged
                    vec=rn(n1s1,1:3)-rn(n2s1,1:3);
                    close_to_n1s1=((L1*L1*(vec*vec'))<mindist2);
                    close_to_n2s1=(((1-L1)*(1-L1)*(vec*vec'))<mindist2);
                    % if collision point is close to one of the existing nodes use that node
                    if close_to_n1s1
                        mergenode2=n1s1;
                    elseif close_to_n2s1
                        mergenode2=n2s1;
                    else
                        spnode=n1s1;
                        splitconnection=linksinconnect(linkid,1);
                        posvel=rn(n1s1,1:lrn3).*(1-L1)+rn(n2s1,1:lrn3).*L1;
                        [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel);
                        mergenode2=length(rn(:,1));
                        newlink=length(links(:,1));
                        links(newlink,6:8)=links(linkid,6:8);
                        fseg=[fseg;zeros(1,6)];
                    end
                    %merge the two nodes
                    disp(sprintf('node %d and node %d are colliding by hinge condition',mergenode2,mergenode1))
                    collisionpoint=findcollisionpoint(mergenode1,mergenode2,rn,connectivity,links);
                    rn(mergenode1,1:lrn2)=[collisionpoint 0 0 0 max(rn(mergenode1,lrn2),rn(mergenode2,lrn2)) ];
                    [rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,mergenode1,mergenode2,MU,NU,a,Ec);
                    if mergednodeid>0
                        for k=1:connectivity(mergednodeid,1)
                            linkid=connectivity(mergednodeid,2*k);
                            fseg(linkid,:)=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),links,appliedstress,linkid);
                            othernode=links(linkid,3-connectivity(mergednodeid,2*k+1));
                            clist=[connectivity(othernode,1) linspace(1,connectivity(othernode,1),connectivity(othernode,1))];
                            [rn(othernode,4:6),fntmp]=feval(mobility,fseg,rn,links,connectivity,othernode,clist);
                        end
                        numbcon=connectivity(mergednodeid,1);
                        conlist=[numbcon linspace(1,numbcon,numbcon)];
                        [rn(mergednodeid,4:6),fntmp]=feval(mobility,fseg,rn,links,connectivity,mergednodeid,conlist);
                    end 
                    %there has been a connectivity change in node i start the search through node i's connections from the beginning
                    if i>size(rn,1)
                        % this is a rare but possible case.
                        % for this condition to be satisfied the last node was being checked for closed hinge condition and it merged with another node
                        % since this was the last node being checked exit the function
                        return;
                    else
                        j=0;
                        k=connectivity(i,1);
                    end
                end
            end
            k=k+1;
        end
        j=j+1;
    end
    i=i+1;
end



function [dist2,ddist2dt,L1,L2]=mindistcalc(x0vx0,x1vx1,y0vy0,y1vy1);
% this function finds the minimum distance bewtween two line segments
% seg1=x0->x1 seg2=y0->y1
% dist2 = square of the minimum distance between the two points
% L1 = normalize position on seg1 that is closest to seg2
% L2 = normalized position on seg2 that is closest to seg1
% ddist2dt = time rate of change of the distance between L1 and L2
x0=x0vx0(1:3);
x1=x1vx1(1:3);
y0=y0vy0(1:3);
y1=y1vy1(1:3);
if length(x0vx0)==6
    vx0=x0vx0(4:6);
    vx1=x1vx1(4:6);
    vy0=y0vy0(4:6);
    vy1=y1vy1(4:6);
else
    vx1=zeros(1,3);
    vx0=zeros(1,3);
    vy1=zeros(1,3);
    vy0=zeros(1,3);
end

seg1=x1-x0;
seg2=y1-y0;
vseg1=vx1-vx0;
vseg2=vy1-vy0;

A=seg1*seg1';
B=2*seg1*(x0'-y0');
C=2*seg1*seg2';
D=2*seg2*(y0'-x0');
E=seg2*seg2';
F=x0*x0'+y0*y0';
G=C*C-4*A*E;
eps=1e-12;
if A<eps % seg1 is just a point
    L1=0;
    if E<eps
        L2=0;
    else
        L2=-0.5*D/E;
    end
elseif E<eps % seg2 is just a point
    L2=0;
    if A<eps
        L1=0;
    else
        L1=-0.5*B/A;
    end
elseif abs(G)<eps % lines are parallel
    dist2=[(y0-x0)*(y0-x0)' (y1-x0)*(y1-x0)' (y0-x1)*(y0-x1)' (y1-x1)*(y1-x1)'];
    [mindist2,pos]=min(dist2);
    L1=floor(pos/2);
    L2=mod(pos-1,2);
else
    L2=(2*A*D+B*C)/G;
    L1=0.5*(C*L2-B)/A;
end

% now check to make sure that L2 and L1 are betwen 0 and 1
L1=min(max([L1,0]),1);
L2=min(max([L2,0]),1);

% now calculate the distance^2 and the time rate of change of the distance between the points at L1 and L2
dist2=(x0+seg1.*L1-y0-seg2.*L2)*(x0+seg1.*L1-y0-seg2.*L2)';
ddist2dt=2*((vx0+vseg1.*L1-vy0-vseg2.*L2)*(x0+seg1.*L1-y0-seg2.*L2)');


 function collisionpoint=findcollisionpoint(mergenode1,mergenode2,rn,connectivity,links);
 % this subroutine finds the collision point of two nodes given that there are strict glide plane constraints
 eps=1e-12;
 newplanecondition=0.875;
 p1=rn(mergenode1,1:3);
 p2=rn(mergenode2,1:3);
 Nmat=zeros(3,3);
 Nsize=0;
 vector=zeros(3,1);
 s=size(rn,2);
 if rn(mergenode1,s)==7
     collisionpoint=rn(mergenode1,1:3);
     return;
 elseif rn(mergenode2,s)==7
     collisionpoint=rn(mergenode2,1:3);
     return;
 end
 
 for i=1:connectivity(mergenode1,1)
     if Nsize<3
         linkid=connectivity(mergenode1,2*i);
         connode=links(connectivity(mergenode1,2*i),3-connectivity(mergenode1,2*i+1));
         rt=rn(mergenode1,1:3)-rn(connode,1:3);                                                              
         L=norm(rt);
         linedir=rt./L;
         n1=cross(linedir,links(linkid,3:5));
         n2=cross(linedir,rn(mergenode1,4:6));
         
         if n1*n1'>eps
             plane=n1./norm(n1);
         elseif n2*n2'>eps
             plane=n2./norm(n2);
         end
         
         if ((n1*n1'>eps)|(n2*n2'>eps))
            if Nsize==0
                conditionismet = 1;
            elseif Nsize==1
                conditionismet = ((Nmat(1,:)*plane')^2 < newplanecondition*newplanecondition);
            else 
                detN=det([Nmat(1:2,:);plane]);
                conditionismet = detN*detN > (1-newplanecondition)^4;
            end
            if conditionismet
                Nsize=Nsize+1;
                Nmat(Nsize,:)=plane;
                vector(Nsize)=plane*p1';
            end
        end
     end
 end
 
 for i=1:connectivity(mergenode2,1)
     if Nsize<3
         linkid=connectivity(mergenode2,2*i);
         connode=links(connectivity(mergenode2,2*i),3-connectivity(mergenode2,2*i+1));
         rt=rn(mergenode2,1:3)-rn(connode,1:3);                                                              
         L=norm(rt);
         linedir=rt./L;
         n1=cross(linedir,links(linkid,3:5));
         n2=cross(linedir,rn(mergenode2,4:6));
         
         if n1*n1'>eps
             plane=n1./norm(n1);
         elseif n2*n2'>eps
             plane=n2./norm(n2);
         end
         if ((n1*n1'>eps)|(n2*n2'>eps))
            if Nsize==1
                conditionismet = ((Nmat(1,:)*plane')^2 < newplanecondition*newplanecondition);
            else 
                detN=det([Nmat(1:2,:);plane]);
                conditionismet = detN*detN > (1-newplanecondition)^4;
            end
         
            if conditionismet
                Nsize=Nsize+1;
                Nmat(Nsize,:)=plane;
                vector(Nsize)=plane*p2';
            end
         end         
     end
 end
 
 Matrix=[eye(3) Nmat(1:Nsize,:)';Nmat(1:Nsize,:) zeros(Nsize,Nsize)];
 V=[(rn(mergenode1,1:3)'+rn(mergenode2,1:3)')./2; vector(1:Nsize)];
 res=Matrix\V;
 collisionpoint=res(1:3)';