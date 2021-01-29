function [connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
%the connectivity is basically the inverse of the links array
% for each node it the connectivity tells how many links that node is part of 
% and then lists the link indices of the system
% connectivity(i,1)= number of links that node i is part of 
% connectivity(i,2:7)= the list of connectivity(i,1) links that node i is part of
rnlength=length(rn(:,4));
connectivity=zeros(rnlength,1+2*maxconnections);
linkslength=length(links(:,1));
linksinconnect=zeros(linkslength,2);
for i=1:linkslength
    if links(i,1)~=0
        a=links(i,1);
        b=links(i,2);
        connectivity(a,1)=connectivity(a,1)+1;
        connectivity(b,1)=connectivity(b,1)+1;
        connectivity(a,2*connectivity(a,1):2*connectivity(a,1)+1)=[i 1];
        connectivity(b,2*connectivity(b,1):2*connectivity(b,1)+1)=[i 2];
        linksinconnect(i,1)=connectivity(a,1);
        linksinconnect(i,2)=connectivity(b,1);
    end
end