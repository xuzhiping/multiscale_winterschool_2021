function Abaqus2PhasefieldUEL(inputPath,MatProp)
narginLoc=nargin;
if narginLoc<2
    if narginLoc==0
        [FileName,PathName,FilterIndex] = uigetfile('*.inp');
        inputPath=[PathName FileName];
    else
    end
    if isa(inputPath,'char')~=1
        error('Provide a valid input file.')
    end
    nMatini=input('Number of materials? - ');
        for i=1:nMatini
            MatProp(i,1)=input(['Young''s modulus for material (' num2str(i) ')? - ']);
            MatProp(i,2)=input(['Poisson''s ratio for material (' num2str(i) ')? - ']);
            MatProp(i,3)=input(['Element thickness for material (' num2str(i) ')? - ']);
            MatProp(i,4)=input(['Stability parameter for material (' num2str(i) ')? - ']);
            MatProp(i,5)=input(['Length scale parameter for material (' num2str(i) ')? - ']);
            MatProp(i,6)=input(['Fracture surface energy release rate for material (' num2str(i) ')? - ']);
        end
end

if length(MatProp(1,:))<6
    error('Not enought material properties are given.')
end
if length(MatProp(1,:))>6
    error('Too much material properties are given.')
end
    
if isempty(inputPath)==1
    [FileName,PathName,FilterIndex] = uigetfile('*.inp');
    inputPath=[PathName FileName];
end

% ---------------- Uploaing material properties -------------------
Ymod=MatProp(:,1);
nu=MatProp(:,2);
thck=MatProp(:,3);
park=MatProp(:,4);

lc=MatProp(:,5);
gc=MatProp(:,6);

% ---------------- Initializing output file -------------------
inputName=inputPath;
output=[inputPath(1:end-4) '_UEL.inp'];

fid=fopen(inputName,'rt');
fout1 = fopen(output, 'wt');
a=fgets(fid);

cnt=0;

% ---------------- Sart reading input file -------------------
cntE=1;
while(ischar(a))

% ------ Repeating nodes ------
while(ischar(a))
    if strfind(a,'*Element')~=0
        break
    end
    fprintf(fout1,a);    
    a=fgets(fid);
end

a=fgets(fid);

% ------ Defintion of UEL types ------
if length(str2num(a))==4
   % triangular elements 
    fprintf(fout1,['*********************** TRIANGULAR ****************************\n']);
    fprintf(fout1,['*User element, nodes=3, type=U1, properties=3, coordinates=2, VARIABLES=6\n']);
    fprintf(fout1,['3\n']);
    fprintf(fout1,['***************************************************************\n']);
    fprintf(fout1,['*User element, nodes=3, type=U3, properties=4, coordinates=2, VARIABLES=42\n']);
    fprintf(fout1,['1,2\n']);
    NnodeE(cntE)=3;
else
   % rectangular elements 
    fprintf(fout1,['*********************** RECTENGULAR ***************************\n']);
    fprintf(fout1,['*User element, nodes=4, type=U2, properties=3, coordinates=2, VARIABLES=8\n']);
    fprintf(fout1,['3\n']);
    fprintf(fout1,['***************************************************************\n']);
    fprintf(fout1,['*User element, nodes=4, type=U4, properties=4, coordinates=2, VARIABLES=56\n']);
    fprintf(fout1,['1,2\n']);
    NnodeE(cntE)=4;
end

% ------ Uploading original mesh ------
cnt=0;
while isempty(strfind(a,'*Nset,')) && isempty(strfind(a,'*Element,'))
    cnt=cnt+1;
    ElMat{cntE}(cnt,:)=str2num(a);
    a=fgets(fid);
end
nElem(cntE)=length(ElMat{cntE}(:,1));

fprintf(fout1,['***************************************************************\n']);

if isempty(strfind(a,'*Element,'))==0
    cntE=cntE+1;
else
    break
end
end

% ------ All elements ------
nElemAll=sum(nElem);

% ------ Replicating UELs ------
for j=1:cntE
    if NnodeE(j)==3
        ej=1;
    elseif NnodeE(j)==4
        ej=2;
    end
    fprintf(fout1,['***************************************************************\n']);
    fprintf(fout1,['*Element, type=U' num2str(ej) '\n']);
    for i=1:nElem(j)
        fprintf(fout1,[num2str(ElMat{j}(i,1)) ', ' num2str(ElMat{j}(i,2:end-1),'%d, ') ' ' num2str(ElMat{j}(i,end)) '\n']);
    end
    fprintf(fout1,['***************************************************************\n']);
    fprintf(fout1,['*Element, type=U' num2str(ej+2) '\n']);
    for i=1:nElem(j)
        fprintf(fout1,[num2str(ElMat{j}(i,1)+nElemAll) ', ' num2str(ElMat{j}(i,2:end-1),'%d, ') ' ' num2str(ElMat{j}(i,end)) '\n']);
    end
end

% ------ Sets and material properties ------
fprintf(fout1,['***************************************************************\n']);
fprintf(fout1,['********************** ASSIGNING MATERIAL PROP ****************\n']);

cnt2=1;
while isempty(strfind(a,'** Section'))
   while isempty(strfind(a,'*Elset,'))
       a=fgets(fid);
   end
   bin=strfind(a,'elset=');
   Ename=a(bin+6:end-1);
   if isempty(strfind(Ename,','))==0
      bin=strfind(Ename,',');
      Ename=Ename(1:bin-1);
   end
   ElstName{cnt2}=Ename;
   fprintf(fout1,a);
   ElSeT{1,cnt2}=a;
   a=fgets(fid);
   cnt3=2;
   while isempty(strfind(a,'** Section:')) && isempty(strfind(a,'*Nset,'))
       ElSeT{cnt3,cnt2}=a;
       fprintf(fout1,a);
       a=fgets(fid);
       cnt3=cnt3+1;
   end
   Nline(cnt2)=cnt3-1;
   cnt2=cnt2+1;
end
fprintf(fout1,['***************************************************************\n']);

nMat=length(ElstName);
if length(Ymod)<nMat
    error('Model contains more materials. Provide more material properties!')
end

for k=1:nMat
    fprintf(fout1,['*Uel property, elset=' ElstName{k} '\n']);
    fprintf(fout1,[num2str(lc(k)) ', ' num2str(gc(k)) ', ' num2str(thck(k)) '\n']);
end
fprintf(fout1,['**lc, gc, h\n']);

fprintf(fout1,['***************************************************************\n']);

for k=1:nMat
    b=ElSeT{1,k};
    c=length(ElstName{k});
    fprintf(fout1,[b(1:14) ElstName{k} '_SS' b(14+c+1:end)]);
    if isempty(strfind(b,'generate'))
        for j=2:Nline(k)-1
            b=str2num(ElSeT{j,k});
            fprintf(fout1,'%d, ',b+nElemAll);
            fprintf(fout1,'\n');
        end
        b=str2num(ElSeT{Nline(k),k});
        fprintf(fout1,'%d, ',b(1:end-1)+nElemAll);
        fprintf(fout1,'%d',b(end)+nElemAll); 
        fprintf(fout1,'\n');
        
    else
        b=str2num(ElSeT{2,k});
        fprintf(fout1,[num2str(b(1)+nElemAll) ', ' num2str(b(2)+nElemAll) ', ' num2str(b(3)) '\n']);
    end
    
end

fprintf(fout1,['***************************************************************\n']);
for k=1:nMat
    fprintf(fout1,['*Uel property, elset=' ElstName{k} '_SS\n']);
    fprintf(fout1,[num2str(Ymod(k)) ', ' num2str(nu(k)) ', ' num2str(thck(k)) ', ' num2str(park(k)) '\n']);
end
fprintf(fout1,['**E, nu, h, k\n']);

% ------ Replicating dummy UMAT elements ------
fprintf(fout1,['***************************************************************\n']);

for j=1:cntE
    fprintf(fout1,['*Element, type=CPE' num2str(NnodeE(j)) '\n']);
    for i=1:nElem(j)
        fprintf(fout1,[num2str(ElMat{j}(i,1)+nElemAll*2) ', ' num2str(ElMat{j}(i,2:end-1),'%d, ') ' ' num2str(ElMat{j}(i,end)) '\n']);
    end
    fprintf(fout1,['***************************************************************\n']);
end

while isempty(strfind(a,'material='))
    a=fgets(fid);
end

UMATname=a(strfind(a,'material=')+9:end);

fprintf(fout1,['*Elset, elset=umatelem, generate\n']);
fprintf(fout1,[num2str(nElemAll*2+1) ', ' num2str(nElemAll*3)   ', 1\n']);

fprintf(fout1,['***************************************************************\n']);
fprintf(fout1,['*Solid Section, elset=umatelem, material=' UMATname '1.0\n']);

while isempty(strfind(a,'*End Part'))
    a=fgets(fid);
end

while(ischar(a))
    if strfind(a,'*Instance, name=')~=0
        break
    end
    fprintf(fout1,a);    
    a=fgets(fid);
end

InstanceName=a(strfind(a,'*Instance, name=')+16:strfind(a,', part=')-1);

while(ischar(a))
    if strfind(a,'*End Assembly')~=0
        break
    end
    fprintf(fout1,a);    
    a=fgets(fid);
end

fprintf(fout1,['**\n']);

fprintf(fout1,['*Elset, elset=umatelem, instance=' InstanceName ', generate\n']);
fprintf(fout1,[num2str(nElemAll*2+1) ', ' num2str(nElemAll*3) ', 1\n']);
fprintf(fout1,['**\n*End Assembly\n']);
fprintf(fout1,['***************************************************************\n']);

a=fgets(fid);

while(ischar(a))

    while(ischar(a))
        if strfind(a,'*Output, history')~=0
            break
        end
        fprintf(fout1,a);    
        a=fgets(fid);
    end

    fprintf(fout1,['*element output, elset=umatelem\n']);
    fprintf(fout1,['SDV\n']);
    fprintf(fout1,['**\n']);
    a=fgets(fid);
    fprintf(fout1,a);
    a=fgets(fid);

end
fclose('all');

% --------------- Start modifying the UEL fortran file -------------------

% determining if the OS if Lunix, Windows of Mac
comp = computer;
if strfind(comp,'WIN')~=0
    comptype=1;     % windows
elseif strfind(comp,'LNX')~=0
    comptype=2;     % linux
elseif strfind(comp,'MAC')~=0
    comptype=3;     % mac
else
    comptype=input('Please provide OS type: 1 - Microsoft Windows; 2 - Linux; 3 - Apple Mac: ');     % asking user
    if comptype==1 || comptype==2 || comptype==3
    else
        comptype=input('Wrong type.\n Please provide OS type: 1 - Microsoft Windows; 2 - Linux; 3 - Apple Mac: ');     % asking user
        if comptype==1 || comptype==2 || comptype==3
        else
        error('Wrong type. Try again.')
    end
    end
end

fullp=mfilename('fullpath');
if isempty(fullp)~=1
    input2=[fullp(1:end-22) '\PhaseFieldUEL\SingleCrack.for'];
else
    fullp=matlab.desktop.editor.getActiveFilename;
    input2=[fullp(1:end-24) '\PhaseFieldUEL\SingleCrack.for'];
end

if comptype==1
    output=[inputName(1:end-4) '_UEL_PS.for'];
else
    output=[inputName(1:end-4) '_UEL_PS.f'];
end

fid=fopen(input2,'rt');
fout1 = fopen(output, 'wt');
a=fgets(fid);

while(ischar(a))
    if strfind(a,'N_ELEM=')~=0
        strs=strfind(a,'N_ELEM=');
        ends=strs+6;
        a=[a(1:ends) num2str(nElemAll) a(ends+2:end)];
    end
    fprintf(fout1,a);    
    a=fgets(fid);
end

fclose('all');

return