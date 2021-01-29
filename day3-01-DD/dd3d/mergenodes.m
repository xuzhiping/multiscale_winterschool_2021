function [rn,connectivity,links,linksinconnect,fseg,nodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,nodeid,deadnode,MU,NU,a,Ec)

% the purpose of this subroutine is to merge the connectivity information in nodeid with deadnode
% then remove deadnode and repair the connectivity list and the link list 
% so that the there are no self links and no two nodes are linked more than once
if nodeid==deadnode
    % check to make sure that the nodes to be merged are not the same node
    % if they are the same node an error has been made in the function call
    % just return without performing any nodal manipulations
    return;
end
nodeid_currentconnections=connectivity(nodeid,1);
deadnode_currentconnections=connectivity(deadnode,1);
total_connections=nodeid_currentconnections+deadnode_currentconnections;
% give the deadnode connections to nodeid
connectivity(nodeid,2*(nodeid_currentconnections+1):2*total_connections+1)=connectivity(deadnode,2:2*deadnode_currentconnections+1);
connectivity(nodeid,1)=total_connections;
% replace deadnode with nodeid where it appears in links
% update linksinconnect to reflect the new positions of the links in the connectivity list
for i=1:deadnode_currentconnections
    links(connectivity(deadnode,2*i),connectivity(deadnode,2*i+1))=nodeid;
    linksinconnect(connectivity(deadnode,2*i),connectivity(deadnode,2*i+1))=nodeid_currentconnections+i;
end
% clear the connectivity list of the deadnode to prepare it for removal
connectivity(deadnode,1:2*deadnode_currentconnections+1)=zeros(1,2*deadnode_currentconnections+1);
% remove the deadnode that is no longer participating in the calculation
[rn,connectivity,links]=removedeadnode(rn,connectivity,links,deadnode);
%if nodeid was the last node in the node list update its position due to execution of the removedeadnode function
if nodeid>length(rn(:,4)) 
    nodeid=deadnode;
end
i=1;
% delete the self links
% self links are links that connect the same node to itself
while i<connectivity(nodeid,1)
    linkid=connectivity(nodeid,2*i);
    posi_linkid=connectivity(nodeid,2*i+1);
    posnoti_linkid=3-posi_linkid;
    nodenoti_linkid=links(linkid,posnoti_linkid);
    if nodeid==nodenoti_linkid
        [links,connectivity,linksinconnect,fseg]=removelink(links,connectivity,linksinconnect,fseg,linkid);
    else
        i=i+1; 
    end
end
% two nodes cannot be connected by more than one link
% search through the connectivity of nodeid and delete multiple links to other nodes
i=1;
while i<connectivity(nodeid,1)
    linkid=connectivity(nodeid,2*i);
    posi_linkid=connectivity(nodeid,2*i+1);
    posnoti_linkid=3-posi_linkid;
    nodenoti_linkid=links(linkid,posnoti_linkid);
    j=i+1;
    while j<=connectivity(nodeid,1)
        linkid2=connectivity(nodeid,2*j);
        posi_linkid2=connectivity(nodeid,2*j+1);
        posnoti_linkid2=3-posi_linkid2;
        nodenoti_linkid2=links(linkid2,posnoti_linkid2);
        if nodenoti_linkid==nodenoti_linkid2
            % a double link has been found
            % create the junction by adjusting the burgers vector of linkid and removing linkid2
            rn0=rn(links(linkid,1),1:3);
            rn1=rn(links(linkid,2),1:3);
            rn2=rn0;
            link1=links(linkid,:);
            link2=links(linkid2,:);
            if posi_linkid+posi_linkid2==3
                links(linkid,3:5)=links(linkid,3:5)-links(linkid2,3:5);
            else %posi_linkid+posi_linkid2~=3
                links(linkid,3:5)=links(linkid,3:5)+links(linkid2,3:5);
                link2=[link2(2) link2(1) -link2(3:5) link2(6:8)];
            end
            % must now change the glide plane of the resultant dislocation
            vec1=rn(nodeid,1:3)-rn(nodenoti_linkid,1:3); % line direction
            vec2=rn(nodeid,4:6)+rn(nodenoti_linkid,4:6); % velocity direction
            n1=cross(vec1,links(linkid,3:5)); % normal determined by line direction and burgers vector  not good for screw dislocations
            n2=cross(vec1,vec2);  % normal determined by line direction and velocity vector this option is used if a screw is found
            if n1*n1'>eps
                links(linkid,6:8)=n1./norm(n1);
            elseif n2*n2'>eps
                links(linkid,6:8)=n2./norm(n2);
            end
            % must remove the second link
            [links,connectivity,linksinconnect,fseg]=removelink(links,connectivity,linksinconnect,fseg,linkid2);
            % if linkid was the last link in the links it position must be updated due to the execution of removelink
            if linkid>length(links(:,1)) 
                linkid=linkid2;
            end
            % check to see if the junction formed by the multiple links has a zero burgers vector
            % if zero bugers vector is found remove the link
            if links(linkid,3:5)*links(linkid,3:5)'==0
                [links,connectivity,linksinconnect,fseg]=removelink(links,connectivity,linksinconnect,fseg,linkid);
                % remove the node that was connected to nodeid if it no longer has any connections
                if connectivity(nodenoti_linkid,1)==0 
                    [rn,connectivity,links]=removedeadnode(rn,connectivity,links,nodenoti_linkid);
                    % if nodeid was the last node in the node list then its position has changed due to execution of removedeadnode
                    if nodeid>length(rn(:,4))
                        nodeid=nodenoti_linkid;
                    end
                end
            end
            % there has been a configurational change so check through the current connection position again
            i=i-1;
            j=connectivity(nodeid,1)+1; 
        else %nodenoti_linkid~=nodenoti_linkid2
            j=j+1;
        end
    end %while j<connectivity(nodeid,1)
    i=i+1;
end %i<connectivity(nodeid,1)
% if nodeid no longer is part of the simulation remove it from the node list
if connectivity(nodeid,1)==0 
    [rn,connectivity,links]=removedeadnode(rn,connectivity,links,nodeid);
    nodeid=0;
end
%disp(sprintf('node %d and node %d have successfully merged',nodeid,deadnode))
% nodes have been successfully merged

function [links,connectivity,linksinconnect,fseg]=removelink(links,connectivity,linksinconnect,fseg,linkid);
% this subroutine is called by meshcoarsen
% this subroutine deletes the link information from the connectivity list
% removes the link from the link list and replaces the linkid with the last link
% after executing this subroutine all of the data structures should be clean

if(linkid>length(links(:,1)))
    disp(sprintf('try to remove link %d while total number of links is %d',linkid,length(links(:,1))));
    pause
end
% delete the linkid where it appears in connectivity list
% first appearance
nodeid1=links(linkid,1);
deadconnection1=linksinconnect(linkid,1);
[connectivity,linksinconnect]=removedeadconnection(connectivity,linksinconnect,nodeid1,deadconnection1);
% second appearance
nodeid2=links(linkid,2);
deadconnection2=linksinconnect(linkid,2);
[connectivity,linksinconnect]=removedeadconnection(connectivity,linksinconnect,nodeid2,deadconnection2);
% remove the link that no longer appears in the connectivity list and doesn't connect any nodes anymore
[links,connectivity,linksinconnect,fseg]=removedeadlink(links,connectivity,linksinconnect,fseg,linkid);

function [linksnew,connectivitynew,linksinconnectnew,fsegnew]=removedeadlink(linksnew,connectivitynew,linksinconnectnew,fsegnew,deadlink);
% this subroutine is called by meshcoarsen
% this subroutine replaces the link in linkid with the link in llinks
% repairs the connectivity and then deletes the llinks from links
% llinks should be the last linkid in links

% change the linkid in the connectivity to reflect the replacement
if(linksinconnectnew(deadlink,:)*linksinconnectnew(deadlink,:)'~=0)
    disp(sprintf('this dead link still has connections and should not be removed'));
    pause
end


llinks=length(linksnew(:,1));
if deadlink<llinks
    linksnew(deadlink,:)=linksnew(llinks,:);
    fsegnew(deadlink,:)=fsegnew(llinks,:);
    linksinconnectnew(deadlink,:)=linksinconnectnew(llinks,:);
    connectivitynew(linksnew(deadlink,1),2*linksinconnectnew(deadlink,1))=deadlink;
    connectivitynew(linksnew(deadlink,2),2*linksinconnectnew(deadlink,2))=deadlink;
end
linksnew(llinks,:)=[];
fsegnew(llinks,:)=[];
linksinconnectnew(llinks,:)=[];

function [rnnew,connectivitynew,linksnew]=removedeadnode(rnnew,connectivitynew,linksnew,deadnode);
% this subroutine is called by meshcoarsen and by removenode
% it removes nodes that are no longer part of the simulation
% and cleans up the data structures
% this subroutine replaces the node in i with the node in lrn
% repairs the links and then deletes the lrn from the node list
% lrn should be the last nodeid in rn
if(connectivitynew(deadnode,1)~=0)
    disp(sprintf('this deadnode still has connections and should not be removed'));
    deadnode
    rnnew
    connectivitynew
    pause
end

lrn=length(rnnew(:,1));
if deadnode<lrn
    rnnew(deadnode,:)=rnnew(lrn,:);
    connectivitynew(deadnode,:)=connectivitynew(lrn,:);
    for j=1:connectivitynew(deadnode,1) % change the nodeid in linksnew from lrn to i
        linksnew(connectivitynew(deadnode,2*j),connectivitynew(deadnode,2*j+1))=deadnode;
    end
end
rnnew(lrn,:)=[];
connectivitynew(lrn,:)=[];

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
               
