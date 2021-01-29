function plotnodes(rn,links,plim)
%plot nodes
%only those nodes within [-plim,plim] in x,y,z directions are plotted

plot3(0,0,0); hold on;
LINKMAX=length(links(:,1));
for i=1:LINKMAX,
    n0=links(i,1);
    n1=links(i,2);
    if((n0~=0)&(n1~=0)&(max(max(abs(rn([n0,n1],:))))<=plim))
        %filter out "infinity" lines
        plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3),'r.-');
    end 
end
hold off
axis equal
grid