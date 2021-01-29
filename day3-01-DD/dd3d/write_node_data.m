%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: write_node_data.m
%
% Purpose:  Write dislocation nodal data into ParaDiS format
%
% Note: rn has to be a compact array (no deleted nodes in the middle)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_node_data(fname,rn,links,L)

%settings
version_number = 2;
filesegments_number = 1;

xBoundMin = -L/2;
yBoundMin = -L/2;
zBoundMin = -L/2;
xBoundMax =  L/2;
yBoundMax =  L/2;
zBoundMax =  L/2;

domainID = 0;
Nnodes = length(rn(:,1));

myflag = 0;

%internal data structure
[NMAX,m]=size(rn);
[LINKMAX,m]=size(links);

%number of nodes
Nnodes = sum((rn(:,4)~=-1));

if(Nnodes ~= NMAX)
  disp(sprintf('rn contains deleted nodes! Nnodes = %d NMAX = %d',Nnodes,NMAX));
end

%build link list
list=zeros(Nnodes,100);
for j=1:LINKMAX,
    n0=links(j,1);
    n1=links(j,2);
    if(n0~=0)&(n1~=0)
        list(n0,1)=list(n0,1)+1;
        list(n0,list(n0,1)+1)=j;
        list(n1,1)=list(n1,1)+1;
        list(n1,list(n1,1)+1)=j;
    end
end

fid = fopen(fname,'w');

%disp(sprintf('data file %s created',fname));
  
fprintf (fid, '#\n');
fprintf (fid, '#	ParaDiS nodal data file (by write_loop_data.m) \n');
fprintf (fid, '#\n');
fprintf (fid, ' \n');
fprintf (fid, '#\n');
fprintf (fid, '#	File version number\n');
fprintf (fid, '#\n');
fprintf (fid, ' \n');
fprintf (fid, '%d\n', version_number);
fprintf (fid, ' \n');
fprintf (fid, '#\n');
fprintf (fid, '#	Number of data file segments\n');
fprintf (fid, '#\n');
fprintf (fid, ' \n');
fprintf (fid, '%d\n', filesegments_number);
fprintf (fid, ' \n');
fprintf (fid, '#\n');
fprintf (fid, '#	Minimum coordinate values (x, y, z)\n');
fprintf (fid, '#\n');
fprintf (fid, '\n');
fprintf (fid, '%d %d %d\n', xBoundMin, yBoundMin, zBoundMin);
fprintf (fid, ' \n');
fprintf (fid, '#\n');
fprintf (fid, '#	Maximum coordinate values (x, y, z)\n');
fprintf (fid, '#\n');
fprintf (fid, ' \n');
fprintf (fid, '%d %d %d\n', xBoundMax, yBoundMax, zBoundMax);
fprintf (fid, ' \n');
fprintf (fid, '#\n');
fprintf (fid, '#	Node count\n');
fprintf (fid, '#\n');
fprintf (fid, ' \n');
fprintf (fid, '%d\n', Nnodes);
fprintf (fid, ' \n');
fprintf (fid, '#\n');
fprintf (fid, '#	Domain geometry (x, y, z)\n');
fprintf (fid, '#\n');
fprintf (fid, ' \n');
fprintf (fid, '1   1   1  \n');
fprintf (fid, ' \n');
fprintf (fid, '#\n');
fprintf (fid, '#	Domain boundaries (x, y, z)\n');
fprintf (fid, '#\n');
fprintf (fid, ' \n');
fprintf (fid, '   %d\n', xBoundMin);
fprintf (fid, '       %d\n', yBoundMin);
fprintf (fid, '           %d\n', zBoundMin);
fprintf (fid, '           %d\n', zBoundMax);
fprintf (fid, '       %d\n', yBoundMax);
fprintf (fid, '   %d\n', xBoundMax);
fprintf (fid, '#\n');
fprintf (fid, '#	Primary lines: node_tag, x, y, z, num_arms, constraint\n');
fprintf (fid, '#	Secondary lines: arm_tag, burgx, burgy, burgz, nx, ny, nz\n');
fprintf (fid, '#\n');
fprintf (fid, '#       length in unit of burgMag\n');
fprintf (fid, ' \n');


for i = 1:Nnodes,
    numNbrs = list(i,1);
    myflag = rn(i,end);
    fprintf (fid, '     %d,%d      %22.14e %22.14e %22.14e %d %d\n', ...
             domainID,i-1,rn(i,1),rn(i,2),rn(i,3),numNbrs,myflag);
    for j = 1:numNbrs,
        k=list(i,j+1);
        n0=links(k,1);
        n1=links(k,2);
        bv=links(k,3:5)';
        nv=links(k,6:8)';
        if(n0==i)
            nb=n1;
        else
            nb=n0;
            bv=-bv;
        end
        fprintf (fid, '           %d,%d      %g %g %g\n', ...
                 domainID,nb-1,bv(1),bv(2),bv(3));
        fprintf (fid, '                 %g %g %g\n', ...
                 nv(1),nv(2),nv(3));
     end
end

fclose(fid);

