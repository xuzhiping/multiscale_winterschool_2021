function [rn,links,connectivity,linksinconnect,fseg]=separation(rn,links,connectivity,linksinconnect,fseg,mobility,MU,NU,a,Ec,mindist,appliedstress);
[lrn,lrn2]=size(rn);

% search through connectivity list
for i=1:lrn
    c=connectivity(i,1);
    if (c>3) & (rn(i,lrn2)==0)
        % disp(sprintf('separation: node %d has %d arms',i,c));
        % a high order node has been found
        % calculate the power dissipated by the connected geometry and not splitting
        ft=zeros(1,3);
        for j=1:c
            linkid=connectivity(i,2*j);
            pos=connectivity(i,2*j+1);
            ft=ft+fseg(linkid,3*(pos-1)+1:3*(pos));
        end
        Powermax=1.05*rn(i,4:6)*ft';
        %initialize the splitting mode that will be undertaken
        splittingmode=0;
        
        % build the list of possible splittingmodes of the multinode
        conlist=buildsplitmodelist(c);
        numsplitmodes=size(conlist,1);
        % conlist is a matrix with dimension numsplitting modes by number of connections
        % it consists of a series of ones and twos detailing which connection in the original
        % connectivity list is connected in the temporary dislocation structure that is being considered
        
        % save the current configuration of the node to be considered for splitting
        refposveli=rn(i,1:6);
        refconnecti=connectivity(i,:);
        for j=1:c
            reffsegi(j,:)=fseg(refconnecti(2*j),:);
        end
        reflinksinconnect=linksinconnect;
        refrn=rn;
        refconnectivity=connectivity;
        reffseg=fseg;
        reflinks=links;
        % begin investigating the power dissipation of the separated configurations
        for j=1:numsplitmodes
            s=conlist(j,:);
            for k=1:c
                count=c-k+1;
                if s(count)==1
                    s(count)=count;
                else
                    s(count)=[];
                end
            end
            [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,i,s,refposveli);
            lastnode=size(rn,1);
            nodelist=[i,lastnode]';
            cmax=size(connectivity,2);
            clist=zeros(2,cmax);
            clist(:,1)=[connectivity(i,1);connectivity(lastnode,1)];
            clist(1,2:1+connectivity(i,1))= linspace(1,connectivity(i,1),connectivity(i,1));
            clist(2,2:1+connectivity(lastnode,1))= linspace(1,connectivity(lastnode,1),connectivity(lastnode,1));
            % do an evaluataion to find out what the splitting direction is
            [vntmp,fntmp]=feval(mobility,fseg,rn,links,connectivity,nodelist,clist);
            vd1=vntmp(1,:)*vntmp(1,:)';
            vd2=vntmp(2,:)*vntmp(2,:)';
            rn([i lastnode],4:6)=vntmp;
            if vd1>=vd2
                dir=-vntmp(1,:)./sqrt(vd1);
                rn(i,1:3)=rn(i,1:3)-mindist*dir;
            else
                dir=vntmp(2,:)./sqrt(vd2);
                rn(lastnode,1:3)=rn(lastnode,1:3)+mindist*dir;
            end
            % check to see if a junction was needed to repair the connectivity
            % initialize the segementforces if necessary
            totalconnectivity=clist(1,1)+clist(2,1);
            if totalconnectivity>refconnecti(1)
                fseg=[fseg;zeros(1,6)];
            end
            % calculate the segment forces still connected to node i
            for k=1:connectivity(i,1)
                linkid1=connectivity(i,2*k);
                fseg(linkid1,:)=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),links,appliedstress,linkid1);
            end 
            % calculate the segment forces connected to the new node
            for k=1:connectivity(lastnode,1)
                linkid2=connectivity(lastnode,2*k);
                fseg(linkid2,:)=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),links,appliedstress,linkid2);
            end 
            % evaluate the power dissipated by this splitting configuration
            [vntmp,fntmp]=feval(mobility,fseg,rn,links,connectivity,nodelist,clist);
            vdiff=vntmp(2,:)-vntmp(1,:);
            if vdiff*dir'>0
                Powertest=fntmp(1,:)*vntmp(1,:)'+fntmp(2,:)*vntmp(2,:)';
                if Powermax<Powertest
                    Powermax=Powertest;
                    splittingmode=j;
                    splittingvel=vntmp;
                    splittingfseg=[fseg(linkid1,:);fseg(linkid2,:)];
                    splittingfind=[linkid1,linkid2];
                end
                
            end
            % test is complete
            % put the configuration back to the starting position 
            rn(i,1:6)=refposveli;
            [rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,i,lastnode,MU,NU,a,Ec);
            connectivity(i,:)=refconnecti;
            for k=1:c
                linkid=refconnecti(2*k);
                posi=refconnecti(2*k+1);
                fseg(refconnecti(2*k),:)=reffsegi(k,:);
                linksinconnect(linkid,posi)=k;
            end
        end
        % all of the tests have been concludeded if there is a candidate for splitting perform the operation 
        % and initialize the forces and the velocities of the affected segments and nodes
        if splittingmode~=0
            % a viable splitting mode was found
            s=conlist(splittingmode,:);
            for k=1:c
                count=c-k+1;
                if s(count)==1
                    s(count)=count;
                else
                    s(count)=[];
                end
            end
            
            disp(sprintf('separation: node %d has %d arms',i,c));
            disp(sprintf('separation: the following connections have split off'));
            for j=1:length(s)
                disp(sprintf('separation: connection %d',s(j)));
            end
            % get the positions and velocities of the nodes after the split
            vel1=splittingvel(1,:);
            vel2=splittingvel(2,:);
            vd1=vel1*vel1';
            vd2=vel2*vel2';
            if vd1>vd2
                dir=-vel1./sqrt(vd1);
                posvel1=[rn(i,1:3)-(mindist*(1+eps)).*dir vel1];
                posvel2=[rn(i,1:3) vel2];
            else
                dir=vel2./sqrt(vd2);
                posvel1=[rn(i,1:3) vel1];
                posvel2=[rn(i,1:3)+(mindist*(1+eps)).*dir vel2];
            end
            % perform splitting operation
            rn(i,1:6)=posvel1;
            [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,i,s,posvel2);
            lastnode=size(rn,1);
            % initialized the segment forces
            for k=1:length(splittingfind)
                fseg(splittingfind(k),:)=splittingfseg(k,:);
            end
            % recalculate the velocities of the two new nodes and the nodes they are connected to
            for k=1:connectivity(i,1)
                linkid=connectivity(i,2*k);
                othernode=links(linkid,3-connectivity(i,2*k+1));
                clist=[connectivity(othernode,1) linspace(1,connectivity(othernode,1),connectivity(othernode,1))];
                [rn(othernode,4:6),fntmp]=feval(mobility,fseg,rn,links,connectivity,othernode,clist);
            end
            
            for k=1:connectivity(lastnode,1)
                linkid=connectivity(lastnode,2*k);
                othernode=links(linkid,3-connectivity(lastnode,2*k+1));
                clist=[connectivity(othernode,1) linspace(1,connectivity(othernode,1),connectivity(othernode,1))];
                [rn(othernode,4:6),fntmp]=feval(mobility,fseg,rn,links,connectivity,othernode,clist);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [conlist]=buildsplitmodelist(numconnections);
    c=numconnections;
    maxsplitarms=floor(0.5*c);
    % calculate the number of possible splitting modes
    numsplitmodes=0;
    conlist=[];
    count=1;
    for j=2:maxsplitarms
        subnumber=factorial(c)/factorial(c-j)/factorial(j);
        splitmodes=createsplitmodelist(c,j);
        if (2*j==c)% this is a special case that handles mirror symmetry situations
            subnumber=0.5*subnumber;
        end
        for k=1:subnumber
            conlist(count,:)=zeros(1,c);
            for i=1:j
                conlist(count,splitmodes(k,i))=1;
            end
            count=count+1;
        end
        numsplitmodes=numsplitmodes+subnumber;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

function [splitmodes]=createsplitmodelist(nodeconnectivity,numsplitarms)
%initialize the size of splitmodes and fill with zeros
    splitmodes=zeros( factorial(nodeconnectivity) / factorial(nodeconnectivity-numsplitarms) / factorial(numsplitarms) , numsplitarms );
    if (numsplitarms==1)
        % if numsplitarms is 1 then the answer is simply a count from one to nodeconnectivity
        splitmodes=linspace(1,nodeconnectivity,nodeconnectivity)';
    else
        % if numsplitarms is greater than one go into a recursive algorithm to build up the system
        lastone=0;
        for i=1:nodeconnectivity-numsplitarms+1
            % initialize the subslitmodes matrix b
            addlength= factorial(nodeconnectivity-i) / factorial(nodeconnectivity-i-numsplitarms+1) / factorial(numsplitarms-1);
            b=zeros(addlength,numsplitarms-1);
            % calculate the subsplitmodes with a call to splitmodes
            b=createsplitmodelist(nodeconnectivity-i,numsplitarms-1);
            % add the subspitmode list b to the growing splitmode list
            splitmodes(lastone+1:lastone+addlength,1)=i.*ones(addlength,1);
            splitmodes(lastone+1:lastone+addlength,2:numsplitarms)=i+b( 1:addlength , 1:numsplitarms-1 );
            % update the current position of the last valid line of the splitmode system
            lastone=lastone+addlength;
        end       
    end


%total number of splitting options
%4 arm multinode
%may split into 3 two arm splitting modes      0.5 * n! / ( (n-2)! * 2! )       3 total
%5 arm mulitnode
%may split into 10 two arm splitting modes           n! / ( (n-2)! * 2! )      10 total
%6 arm mulitinode
%may split into 15 two arm splitting modes or        n! / ( (n-2)! * 2! )
%may split into 10 three arm spliting modes    0.5 * n! / ( (n-3)! * 3! )      25 total
%7 arm multinode
%may split into 21 two arm splitting modes or        n! / ( (n-2)! * 2! )
%may split into 35 thee arm splitting modes          n! / ( (n-3)! * 3! )      46 total
%8 arm multinode
%may split into 28 two arm splitting modes or        n! / ( (n-2)! * 2! )
%may split into 56 three arm splitting modes or      n! / ( (n-3)! * 3! )
%may split into 35 four arm splitting modes    0.5 * n! / ( (n-4)! * 4! )     119 total 
%9 arm multinode
%may split into 36 two arm splitting modes or        n! / ( (n-2)! * 2! )
%may split into 84 three arm splitting modes or      n! / ( (n-3)! * 3! )
%may split into 126 four arm splitting modes         n! / ( (n-4)! * 4! )     246 total
%10 arm multinode
%may split into 45 two arm splitting modes or        n! / ( (n-2)! * 2! )
%may split into 120 three arm splitting modes or     n! / ( (n-3)! * 3! )
%may split into 210 four arm splitting modes or      n! / ( (n-4)! * 4! )
%may split into 126 five arm splitting modes   0.5 * n! / ( (n-5)! * 5! )     501 total
% the division by 2 for some of the cases is done outside this function call when appropriate.    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


