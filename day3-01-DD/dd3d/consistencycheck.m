function c=consistencycheck(rn,links,connectivity,linksinconnect)
%self consistency check for topology
%
%items to check:
% 1. the length of the connectivity matrix should be the same as that of links
%    (specific to DDLab implementation)
% 2. connectivity matrix: each row correspond to a node, 
%    only the first 2*numNbrs+1 entries in each row can be non-zero
%    where numNbrs is the number of neighbors for this node
%    (specific to DDLab implementation)
% 3. Burgers vector conservation at every node (if flag != 7)
%    (general)
% 4. a node cannot be connected to another node twice
%    (general)
% 5. a link cannot have zero Burgers vector
%    (general)
% 6. a node cannot be connected to the same link twice in the connectivity matrix
%    (specific to DDLab implementation)
% 7. consistency between "connectivity" matrix and "links" matrix
%    if node i has a link j, then link j must have node i as one of its end nodes
%    (specific to DDLab implementation)
% 8. check the consistency between the "connectivity" matrix and the "linkinconnect" matrix
%    (specific to DDLab implementation)

%item 1
% the length of the connectivity matrix should be the same as that of links
if(max(max(connectivity))~=length(links(:,1)))
    disp(sprintf('max connectivity = %d length(links) = %d',max(max(connectivity)),length(links(:,1))))
    pause
end

%item 2
% only the first 2*numNbrs+1 entries in each row can be non-zero
nnodes=length(rn(:,4));
for i=1:nnodes,
    numNbrs=connectivity(i,1);
    if(sum(connectivity(i,2*(numNbrs+1):end)~=0))
        disp('inconsistent connectivity list');
        i
        connectivity(i,:)
        pause
    end
end

%item 3
% Burgers vector conservation at every node
for i=1:nnodes
    if rn(i,4)~=7
        totalb=zeros(1,3);
        numNbrs=connectivity(i,1);
        for j=1:numNbrs
        linkid=connectivity(i,2*j);
        posi=connectivity(i,2*j+1);
        totalb=totalb+(3-2*posi).*links(linkid,3:5);  
        end
        if totalb*totalb'~=0
            disp(sprintf('the Burgers vector is not conserved at node %d',i));
            totalb*totalb'
            pause
        end
    end
end

%if(sum(sum(connectivity==length(links(:,1))))~=2)
%    disp(sprintf('arm %d has %d ends',length(links(:,1)),sum(sum(connectivity==length(links(:,1))))))
%    pause
%end

%item 4
% a node cannot be connected to another node twice
for i=1:nnodes,
    % first find all the nodes connected to node i
    numNbrs=connectivity(i,1);
    neighborNodes = ones(1,numNbrs)*(-1);
    for j=1:numNbrs,
        linkid = connectivity(i,2*j);
        neighborNodes(j) = links(linkid, 3-connectivity(i,2*j+1));
    end
    % then find if there are any repeats in the neighbor nodes
    for j=1:numNbrs-1,
        for k=j+1:numNbrs,
            if neighborNodes(j)==neighborNodes(k)
                disp(sprintf('node %d connected with node %d twice',i,neighborNodes(j)));
                pause
            end
        end
    end
end 
    
%item 5
% a link cannot have zero Burgers vector
for j=1:length(links(:,1))
    b = links(j,3:5);
    if (norm(b)==0)
         disp(sprintf('link %d has zero Burgers vector',j));
         pause
    end
end

%item 6
% a node cannot be connected to the same link twice in the connectivity matrix
for i=1:nnodes,
    numNbrs=connectivity(i,1);
    for j=1:numNbrs,
        for k=j+1:numNbrs,
            if(connectivity(i,2*j)==connectivity(i,2*k))
                connectivity(i,:)
                disp(sprintf('node %d connected with link %d twice',i,connectivity(i,2*j)));
                pause
            end
        end
    end
end


%item 7
% consistency between "connectivity" matrix and "links" matrix
for i=1:nnodes,
    numNbrs=connectivity(i,1);
    for j=1:numNbrs,
        n0=links(connectivity(i,2*j),connectivity(i,2*j+1));    
        if(i~=n0)
            disp(sprintf('connectivity(%d,:)',i));
            connectivity(i,:)
            disp(sprintf('links(%d,:)=(%d,%d)',connectivity(i,2*j),links(connectivity(i,2*j),1),links(connectivity(i,2*j),2)));
        end
    end
end

%item 8
% check the consistency between the "connectivity" matrix and the "linkinconnect" matrix
for i=1:length(links(:,1))    
    j=connectivity(links(i,1),2*linksinconnect(i,1));
    k=connectivity(links(i,2),2*linksinconnect(i,2));
    if(i~=j)|(i~=k)
        disp(sprintf('inconsistent link %d %d %d',i,j,k));
        disp(sprintf('links(%d,:)=(%d,%d)',i,links(i,1),links(i,2)));
        disp(sprintf('linksinconnect(%d,:)=(%d,%d)',i,linksinconnect(i,1),linksinconnect(i,2)));
        disp(sprintf('connectivity(%d,:)=(%d %d %d %d %d %d %d %d %d)',links(i,1),...
            connectivity(links(i,1),1),connectivity(links(i,1),2),connectivity(links(i,1),3),...
            connectivity(links(i,1),4),connectivity(links(i,1),5),connectivity(links(i,1),6),...
            connectivity(links(i,1),7),connectivity(links(i,1),8),connectivity(links(i,1),9) ));
        disp(sprintf('connectivity(%d,:)=(%d %d %d %d %d %d %d %d %d)',links(i,2),...
            connectivity(links(i,2),1),connectivity(links(i,2),2),connectivity(links(i,2),3),...
            connectivity(links(i,2),4),connectivity(links(i,2),5),connectivity(links(i,2),6),...
            connectivity(links(i,2),7),connectivity(links(i,2),8),connectivity(links(i,2),9) ));
            
        pause
    end
end
