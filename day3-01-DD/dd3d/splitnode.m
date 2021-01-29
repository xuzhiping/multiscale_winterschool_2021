function [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,nodeid,splitconnections,posvel);
% this function splits up the connectivity of nodeid with a new node that is added to the end of rn
% after the node is added and the connectivity list is split a check is done to determine whether a link from the new node
% to nodeid is needed to conserve burgers vector. if a link is needed one is added
% add newnode
eps=1e-12;
testl=length(splitconnections);
if ((testl==0)|(testl==connectivity(nodeid,:))|(max(splitconnections)>connectivity(nodeid,:)))
    %sanity check%
    disp('node was not split because of a problem with the connections to be split');
    return;
end
% delete split connections from nodeid
s=splitconnections;
c=connectivity(nodeid,:);
for i=1:length(s)
    % because of reordering that is done remove the biggest connection index first
    [maxs,position]=max(s);
    [connectivity,linksinconnect]=removedeadconnection(connectivity,linksinconnect,nodeid,maxs);
    s(position)=0;
end
newnodeid=size(rn,1)+1;
rn(newnodeid,:)=[posvel 0];
% transfer connectivitylist to be split and update links and linksinconnect
connectivity(newnodeid,1)=length(splitconnections);
for i=1:connectivity(newnodeid,1)
    connectivity(newnodeid,2*i:2*i+1)=c(2*splitconnections(i):2*splitconnections(i)+1);
    links(connectivity(newnodeid,2*i),connectivity(newnodeid,2*i+1))=newnodeid;
    linksinconnect(connectivity(newnodeid,2*i),connectivity(newnodeid,2*i+1))=i;
end
% all of the datastructures should be clean at this point
% now we must check to see if a new link must be made bewteen the newnode and the nodeid for conservation of burgers vector
burgv=links(connectivity(newnodeid,2),3:5);
posnew=connectivity(newnodeid,3);
for i=2:connectivity(newnodeid,1)
    if posnew+connectivity(newnodeid,2*i+1)==3
        burgv=burgv-links(connectivity(newnodeid,2*i),3:5);
    else
        burgv=burgv+links(connectivity(newnodeid,2*i),3:5);
    end
end
if burgv*burgv'~=0
    % must create a link to conserve burgers vector
    lastlink=size(links,1)+1;
    links(lastlink,:)=[ 0 0 burgv 0 0 0];
    links(lastlink,posnew)=nodeid;
    links(lastlink,3-posnew)=newnodeid;
    % update connectivity list of nodeid
    lastconnection1=connectivity(nodeid,1)+1;
    connectivity(nodeid,2*lastconnection1:2*lastconnection1+1)=[lastlink posnew];
    connectivity(nodeid,1)=lastconnection1;
    % update connectivity list of newnodeid
    lastconnection2=connectivity(newnodeid,1)+1;
    connectivity(newnodeid,2*lastconnection2:2*lastconnection2+1)=[lastlink 3-posnew];
    connectivity(newnodeid,1)=lastconnection2;
    % update linksinconnect
    linksinconnect=[linksinconnect; 0 0];
    linksinconnect(lastlink,posnew)=lastconnection1;
    linksinconnect(lastlink,3-posnew)=lastconnection2;
    % calculate the plane of the new link
    % this calculation is suspect because it is unclear what role the plane plays in the mobility function
    % and whether its sign matters
    vec1=rn(nodeid,1:3)-rn(newnodeid,1:3); % line direction
    vec2=rn(nodeid,4:6)+rn(newnodeid,4:6); % velocity direction
    n1=cross(vec1,burgv); % normal determined by line direction and burgers vector  not good for screw dislocations
    n2=cross(vec1,vec2);  % normal determined by line direction and velocity vector this option is used if a screw is found
    if n1*n1'>eps
        links(lastlink,6:8)=n1./norm(n1);
    elseif n2*n2'>eps
        links(lastlink,6:8)=n2./norm(n2);
    end
end
%disp(sprintf('node %d has split into node %d and node %d',nodeid,nodeid,newnodeid));
% end of function splitnode

function [connectivity,linksinconnect]=removedeadconnection(connectivity,linksinconnect,nodeid,deadconnection);
%This subroutine deletes an entry in a node's connectivity list and updates the linksinconnet array

lastconnection=connectivity(nodeid,1);
%remove the entry in linksinconnect to show that the connectivity data no longer exits for that link
linksinconnect(connectivity(nodeid,2*deadconnection),connectivity(nodeid,2*deadconnection+1))=0;
if (lastconnection>deadconnection)
   %replace link in the connectivitylist with the lastlink in the connectivity list
   connectivity(nodeid,2*deadconnection:2*deadconnection+1)=connectivity(nodeid,2*lastconnection:2*lastconnection+1);
   % update linksinconnect to reflect the change in position of the lastconnection
   linksinconnect(connectivity(nodeid,2*deadconnection),connectivity(nodeid,2*deadconnection+1))=deadconnection;
end
connectivity(nodeid,2*lastconnection:2*lastconnection+1)=[0 0];
connectivity(nodeid,1)=lastconnection-1;
               