function [fseg]=segforcevec(MU,NU,a,Ec,rn,links,sigext,linkid);
%compute nodal driving force of dislocation network by stress formula
%(vectorized version)
%rn: nodal position
%links: dislocation segments (connectivity, Burgers vector)
%sigext: external stress (applied)

%segment contribution to node force

%NMAX: max number of nodes[NMAX,m]=size(rn);
[NMAX,m]=size(rn);
if(m~=4)
    disp('rn should have 4 columns!');
    return;
end
[LINKMAX,m]=size(links);

%construct segment list
segments=constructsegmentlist(rn,links);
if linkid==0
    %PK force due to applied stress
    %t0=clock;
    fpk=pkforcevec(sigext,segments);
    %t=etime(clock,t0); disp(sprintf('pkforcevec time = %6.2f seconds\n',t));

    %self force due to self stress
    %t0=clock;
    [fs0,fs1]=selfforcevec(MU,NU,a,Ec,segments);
    %t=etime(clock,t0); disp(sprintf('selfforcevec time = %6.2f seconds\n',t));

    %remote force due to remote stress
    %t0=clock;
    %[fr0,fr1]=remoteforcevec(MU,NU,a,segments,0);
    [fr0,fr1]=remoteforcevec_mex(MU,NU,a,segments,0);
    %t=etime(clock,t0); disp(sprintf('remoteforcevec time = %6.2f seconds\n',t));
    %add force contributions together
    fseg=[fpk, fpk]*0.5+[fs0, fs1]+[fr0, fr1];
else
    %PK force due to applied stress
    %t0=clock;
    fpk=pkforcevec(sigext,segments(linkid,:));
    %t=etime(clock,t0); disp(sprintf('pkforcevec time = %6.2f seconds\n',t));

    %self force due to self stress
    %t0=clock;
    [fs0,fs1]=selfforcevec(MU,NU,a,Ec,segments(linkid,:));
    %t=etime(clock,t0); disp(sprintf('selfforcevec time = %6.2f seconds\n',t));

    %remote force due to remote stress
    %t0=clock;
    %[fr0,fr1]=remoteforcevec(MU,NU,a,segments,linkid);
    [fr0,fr1]=remoteforcevec_mex(MU,NU,a,segments,linkid);
    %t=etime(clock,t0); disp(sprintf('remoteforcevec time = %6.2f seconds\n',t));
    %add force contributions together
    fseg=[fpk, fpk]*0.5+[fs0, fs1]+[fr0, fr1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segments=constructsegmentlist(rn,links)

[LINKMAX,m]=size(links);

segments=zeros(LINKMAX,14);
nseg=0;
for i=1:LINKMAX,
    n0=links(i,1);
    n1=links(i,2);
    if((n0~=0)&(n1~=0))
        nseg=nseg+1;
        segments(nseg,:)=[links(i,1:5),rn(n0,1:3),rn(n1,1:3),links(i,6:8)];
    end
end
segments=segments(1:nseg,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f=pkforcevec(sigext,segments)
%nodal force on dislocation segments due to applied stress sigext
%(vectorized version)
%format of segments array:
% (n0,n1,bx,by,bz,x0,y0,z0,x1,y1,z1,nx,ny,nz)
[nseg,m]=size(segments);
f=zeros(nseg,3);

b=segments(:,3:5);
sigb=sigext*b';
r0=segments(:,6:8);
r1=segments(:,9:11);
r01=r1-r0;

for i=1:nseg,
    f(i,:)=cross(sigb(:,i)',r01(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f1,f2]=selfforcevec(mu,nu,a,Ec,segments)
%self force (due to self stress)
%(vectorized version)
%format of segments array:
% (n0,n1,bx,by,bz,x0,y0,z0,x1,y1,z1)


    Diff=segments(:,9:11)-segments(:,6:8);
 
    L=sqrt(sum(Diff.*Diff,2));
    Linv=1./L;
    
    La=sqrt(L.*L+a*a);
    Lainv=1./La;
    
    t=Diff.*[Linv Linv Linv];
    
    omninv=1/(1-nu);
    
    bs=sum(segments(:,3:5).*t,2);
    bs2=bs.*bs;
    bev=segments(:,3:5)-[bs bs bs].*t;
    be2=sum(bev.*bev,2);
    
    % Elastic Self Interaction Force - Torsional Contribution
    S=(0.25*mu/pi).*bs.*((nu*omninv).*( log((La+L)./a)- 2.*(La-a).*Linv)- 0.5.*(La-a).*(La-a).*Linv.*Lainv);
    % Core Self Interaction Force - Torsional Contribution
    Score=2.*nu*omninv*Ec.*bs;
    Stot=S+Score;
    f2=[Stot Stot Stot].*bev;
    
    % Core Self Interaction Force - Longitudinal Component
    LTcore=(bs2 + be2.*omninv).*Ec;
    f2=f2-[LTcore LTcore LTcore].*t;
    f1=-f2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f0,f1]=remoteforcevec_mex(MU,NU,a,segments,linkid)
%nodal force on dislocation segment 01 due to another segment 23
%format of segments array:
% (n0,n1,bx,by,bz,x0,y0,z0,x1,y1,z1,,nx,ny,nz)
lseg=size(segments,1);
if linkid==0
    f0=zeros(lseg,3);
    f1=zeros(lseg,3);
    k=1;
    for i=1:lseg
        b1=segments(i,3:5);
        x1=segments(i,6:8);
        x2=segments(i,9:11);
        for j=(i+1):lseg
            b2=segments(j,3:5);
            x3=segments(j,6:8);
            x4=segments(j,9:11);
            p1x=x1(1);
            p1y=x1(2);
            p1z=x1(3);
            p2x=x2(1);
            p2y=x2(2);
            p2z=x2(3);
            bpx=b1(1);
            bpy=b1(2);
            bpz=b1(3);
            p3x=x3(1);
            p3y=x3(2);
            p3z=x3(3);
            p4x=x4(1);
            p4y=x4(2);
            p4z=x4(3);
            bx=b2(1);
            by=b2(2);
            bz=b2(3);
            seg12Local=1;
            seg34Local=1;
            [f1x,f1y,f1z,...
                    f2x,f2y,f2z,...
                    f3x,f3y,f3z,...
                    f4x,f4y,f4z]=SegSegForces(p1x,p1y,p1z,...
                p2x, p2y, p2z,...
                p3x, p3y, p3z,...
                p4x, p4y, p4z,...
                bpx, bpy, bpz, ...
                    bx, by, bz,...
                a, MU, NU, seg12Local, seg34Local);
            f0(i,1)=f0(i,1)+f1x;
            f0(i,2)=f0(i,2)+f1y;
            f0(i,3)=f0(i,3)+f1z;
            
            f0(j,1)=f0(j,1)+f3x;
            f0(j,2)=f0(j,2)+f3y;
            f0(j,3)=f0(j,3)+f3z;
            
            f1(i,1)=f1(i,1)+f2x;
            f1(i,2)=f1(i,2)+f2y;
            f1(i,3)=f1(i,3)+f2z;
            
            f1(j,1)=f1(j,1)+f4x;
            f1(j,2)=f1(j,2)+f4y;
            f1(j,3)=f1(j,3)+f4z;
            
        end
    end
else
    tmp=segments(1,:);
    segments(1,:)=segments(linkid,:);
    segments(linkid,:)=tmp;
    f0=zeros(1,3);
    f1=zeros(1,3);
    k=1;
    b1=segments(1,3:5);
    x1=segments(1,6:8);
    x2=segments(1,9:11);
    for j=2:lseg
        b2=segments(j,3:5);
        x3=segments(j,6:8);
        x4=segments(j,9:11);
        p1x=x1(1);
            p1y=x1(2);
            p1z=x1(3);
            p2x=x2(1);
            p2y=x2(2);
            p2z=x2(3);
            bpx=b1(1);
            bpy=b1(2);
            bpz=b1(3);
            p3x=x3(1);
            p3y=x3(2);
            p3z=x3(3);
            p4x=x4(1);
            p4y=x4(2);
            p4z=x4(3);
            bx=b2(1);
            by=b2(2);
            bz=b2(3);
            seg12Local=1;
            seg34Local=0;
            [f1x,f1y,f1z,...
                    f2x,f2y,f2z,...
                    f3x,f3y,f3z,...
                    f4x,f4y,f4z]=SegSegForces(p1x,p1y,p1z,...
                p2x, p2y, p2z,...
                p3x, p3y, p3z,...
                p4x, p4y, p4z,...
                bpx, bpy, bpz, ...
                    bx, by, bz,...
                a, MU, NU, seg12Local, seg34Local);
            f0(1,1)=f0(1,1)+f1x;
            f0(1,2)=f0(1,2)+f1y;
            f0(1,3)=f0(1,3)+f1z;
            f1(1,1)=f1(1,1)+f2x;
            f1(1,2)=f1(1,2)+f2y;
            f1(1,3)=f1(1,3)+f2z;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f0,f1]=remoteforcevec(MU,NU,a,segments,linkid)
%nodal force on dislocation segment 01 due to another segment 23
%format of segments array:
% (n0,n1,bx,by,bz,x0,y0,z0,x1,y1,z1,,nx,ny,nz)
lseg=size(segments,1);
if linkid==0
    f0=zeros(lseg,3);
    f1=zeros(lseg,3);
    k=1;
    for i=1:lseg
        b1=segments(i,3:5);
        x1=segments(i,6:8);
        x2=segments(i,9:11);
        for j=(i+1):lseg
            b2=segments(j,3:5);
            x3=segments(j,6:8);
            x4=segments(j,9:11);
            pairindex(k,:)=[i,j];
            pairs(k,:)=[x1,x2,x3,x4,b1,b2];
            k=k+1;
        end
    end

    np=size(pairs,1);
    [fa,fb,fc,fd]=RemoteNodeForce(pairs(:,1:3),pairs(:,4:6),pairs(:,7:9),pairs(:,10:12),pairs(:,13:15),pairs(:,16:18),a,MU,NU);

    for k=1:np
        i=pairindex(k,1);
        j=pairindex(k,2);
        f0(i,:)=f0(i,:)+fa(k,:);
        f1(i,:)=f1(i,:)+fb(k,:);
        f0(j,:)=f0(j,:)+fc(k,:);
        f1(j,:)=f1(j,:)+fd(k,:);
    end
else
    tmp=segments(1,:);
    segments(1,:)=segments(linkid,:);
    segments(linkid,:)=tmp;
    f0=zeros(1,3);
    f1=zeros(1,3);
    k=1;
    b1=segments(1,3:5);
    x1=segments(1,6:8);
    x2=segments(1,9:11);
    for j=2:lseg
        b2=segments(j,3:5);
        x3=segments(j,6:8);
        x4=segments(j,9:11);
        pairs(k,:)=[x1,x2,x3,x4,b1,b2];
        k=k+1;
    end
    [fa,fb,fc,fd]=RemoteNodeForce(pairs(:,1:3),pairs(:,4:6),pairs(:,7:9),pairs(:,10:12),pairs(:,13:15),pairs(:,16:18),a,MU,NU);
    f0=sum(fa,1);
    f1=sum(fb,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f1,f2,f3,f4]=RemoteNodeForce(x1,x2,x3,x4,bp,b,a,mu,nu);
% this calculates the forces between dislocation nodes analytically
%inputs: endpoints of first dislocation segment starting at x1 ending at x2 with burgers vector bp
%        endpoints of second dislocation segment starting at x3 ending at x4 with burgers vector b
%        core paramter a
%        shear modulus mu
%        poisson ration nu
%
%outputs: f1,f2,f3,f4 is the force on nodes located at x1, x2, x3, x4 respectively
                
    f1=[]; 
    f2=[];
    f4=[];
    f3=[];
    eps=1e-6;
    
    Diff=x4-x3;
    oneoverL=1./sqrt(sum(Diff.*Diff,2));
    t=Diff.*[oneoverL oneoverL oneoverL];
    Diff=x2-x1;
    oneoverLp=1./sqrt(sum(Diff.*Diff,2));
    tp=Diff.*[oneoverLp oneoverLp oneoverLp];

    c=sum(t.*tp,2);
    c2=c.*c;
    onemc2=1-c2;  
    cL=size(c,1);
    k=1;
    spindex=[];
    for i=1:cL
        index=cL+1-i;
        if onemc2(index)<eps
            spindex(k)=index;
            x1sp(k,:)=x1(index,:);
            x2sp(k,:)=x2(index,:);
            x3sp(k,:)=x3(index,:);
            x4sp(k,:)=x4(index,:);
            bsp(k,:)=b(index,:);
            bpsp(k,:)=bp(index,:);
            x1(index,:)=[];
            x2(index,:)=[];
            x3(index,:)=[];
            x4(index,:)=[];
            bp(index,:)=[];
            b(index,:)=[];
            t(index,:)=[];
            tp(index,:)=[];
            oneoverL(index,:)=[];
            oneoverLp(index,:)=[];
            c(index,:)=[];
            c2(index,:)=[];
            onemc2(index,:)=[];
            k=k+1;
        end
    end
    sL=length(spindex);
    if (cL-sL) > 0
        txtp=[t(:,2).*tp(:,3)-t(:,3).*tp(:,2) , t(:,3).*tp(:,1)-t(:,1).*tp(:,3) , t(:,1).*tp(:,2)-t(:,2).*tp(:,1)];
        onemc2inv=1./onemc2;
        R=[ x3-x1 , x4-x2 ];
        d=sum(R(:,1:3).*txtp,2).*onemc2inv;
        temp1=[sum(R(:,1:3).*t,2) sum(R(:,4:6).*t,2)];
        temp2=[sum(R(:,1:3).*tp,2) sum(R(:,4:6).*tp,2)];
        y=(temp1-[ c c ].*temp2).*[ onemc2inv onemc2inv ];
        z=(temp2-[ c c ].*temp1).*[ onemc2inv onemc2inv ];
        
        yin=[y(:,1) y(:,1) y(:,2) y(:,2)];
        zin=[z(:,1) z(:,2) z(:,1) z(:,2)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  this section calculates the formulae from the integral expressions
        a2=a*a;
        a2_d2 = a2+d.*d.*onemc2;
        y2    = yin.*yin;
        z2    = zin.*zin;
        Ra    = sqrt( [a2_d2 a2_d2 a2_d2 a2_d2] + y2 + z2 + 2.*yin.*zin.*[c c c c] );        
        Rainv = 1./Ra;
        
        Ra_Rdot_tp = Ra+zin+yin.*[c c c c];       
        Ra_Rdot_t  = Ra+yin+zin.*[c c c c];       
         
         log_Ra_Rdot_tp =     log(Ra_Rdot_tp);
        ylog_Ra_Rdot_tp = yin.*log_Ra_Rdot_tp;
        
         log_Ra_Rdot_t =     log(Ra_Rdot_t);
        zlog_Ra_Rdot_t = zin.*log_Ra_Rdot_t;
        
          Ra2_R_tpinv = Rainv./Ra_Rdot_tp;
         yRa2_R_tpinv = yin.*  Ra2_R_tpinv;
        y2Ra2_R_tpinv = yin.* yRa2_R_tpinv;
        
          Ra2_R_tinv = Rainv./Ra_Rdot_t;
         zRa2_R_tinv = zin.* Ra2_R_tinv;
        z2Ra2_R_tinv = zin.*zRa2_R_tinv;

        
        denom=1./sqrt(onemc2.*a2_d2);  
        cdenom=(1+c).*denom;
        
        f_003=-2.*[denom denom denom denom].*atan((Ra+yin+zin).*[cdenom cdenom cdenom cdenom]);
        
        adf_003=[a2_d2 a2_d2 a2_d2 a2_d2].*f_003;
        commonf223=( [c c c c].*Ra - adf_003 ).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        
        f_103=( [c c c c].*log_Ra_Rdot_t  - log_Ra_Rdot_tp ).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        f_013=( [c c c c].*log_Ra_Rdot_tp - log_Ra_Rdot_t  ).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        f_113=( [c c c c].*adf_003 - Ra ).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        f_203= zlog_Ra_Rdot_t  + commonf223;
        f_023= ylog_Ra_Rdot_tp + commonf223;
        
        commonf225=f_003 - [c c c c].*Rainv;
        commonf025=[c c c c].*yRa2_R_tpinv - Rainv  ;
        ycommonf025=yin.*commonf025;
        commonf205=[c c c c].*zRa2_R_tinv  - Rainv  ;
        zcommonf205=zin.*commonf205;
        commonf305=log_Ra_Rdot_t  -(yin-[c c c c].*zin).*Rainv - [c2 c2 c2 c2].*z2Ra2_R_tinv;
        zcommonf305=zin.*commonf305;
        commonf035=log_Ra_Rdot_tp -(zin-[c c c c].*yin).*Rainv - [c2 c2 c2 c2].*y2Ra2_R_tpinv;
        tf_113=2.*f_113;
        
        f_005=( f_003 - yRa2_R_tpinv - zRa2_R_tinv )./ [a2_d2 a2_d2 a2_d2 a2_d2];
        f_105=( Ra2_R_tpinv - [c c c c].*Ra2_R_tinv  ).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        f_015=( Ra2_R_tinv  - [c c c c].*Ra2_R_tpinv ).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        f_115=( Rainv - [c c c c].*( yRa2_R_tpinv + zRa2_R_tinv + f_003 )).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        f_205=( yRa2_R_tpinv + [c2 c2 c2 c2].*zRa2_R_tinv  + commonf225 ).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        f_025=( zRa2_R_tinv  + [c2 c2 c2 c2].*yRa2_R_tpinv + commonf225 ).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        f_215=( f_013 - ycommonf025 + [c c c c].*(zcommonf205-f_103) ).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        f_125=( f_103 - zcommonf205 + [c c c c].*(ycommonf025 - f_013) ).*[onemc2inv onemc2inv onemc2inv onemc2inv]; 
        f_225=( f_203 - zcommonf305 + [c c c c].*( y2.*commonf025 - tf_113) ).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        f_305=(y2Ra2_R_tpinv + [c c c c].*commonf305 + 2.*f_103).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        f_035=(z2Ra2_R_tinv  + [c c c c].*commonf035 + 2.*f_013).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        f_315=(tf_113 - y2.*commonf025 + [c c c c].*(zcommonf305 - f_203)).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        f_135=(tf_113 - z2.*commonf205 + [c c c c].*(yin.*commonf035 - f_023)).*[onemc2inv onemc2inv onemc2inv onemc2inv];
        mf=[1 -1 -1 1]';
        Fintegrals=[f_003*mf f_103*mf f_013*mf f_113*mf f_203*mf f_023*mf f_005*mf f_105*mf f_015*mf f_115*mf f_205*mf f_025*mf f_215*mf f_125*mf f_225*mf f_305*mf f_035*mf f_315*mf f_135*mf];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this section calculates the dot products and cross prodcucts for the coefficients
        m4p=0.25 * mu / pi;
        m4pd= m4p .* d;
        m8p=0.5 * m4p;
        m8pd=m8p .* d;
        m4pn=m4p / ( 1 - nu );
        m4pnd=m4pn .* d;
        m4pnd2=m4pnd .* d;
        m4pnd3=m4pnd2 .* d;
        a2m4pnd=a2 .* m4pnd;
        a2m8pd=a2 .* m8pd;
        a2m4pn=a2 * m4pn;
        a2m8p=a2 * m8p;
        
        tpxt=-txtp;
        txbp=[t(:,2).*bp(:,3)-t(:,3).*bp(:,2) , t(:,3).*bp(:,1)-t(:,1).*bp(:,3) , t(:,1).*bp(:,2)-t(:,2).*bp(:,1)];
        tpxb=[tp(:,2).*b(:,3)-tp(:,3).*b(:,2) , tp(:,3).*b(:,1)-tp(:,1).*b(:,3) , tp(:,1).*b(:,2)-tp(:,2).*b(:,1)];
        bxt=[b(:,2).*t(:,3)-b(:,3).*t(:,2) , b(:,3).*t(:,1)-b(:,1).*t(:,3) , b(:,1).*t(:,2)-b(:,2).*t(:,1)];
        bpxtp=[bp(:,2).*tp(:,3)-bp(:,3).*tp(:,2) , bp(:,3).*tp(:,1)-bp(:,1).*tp(:,3) , bp(:,1).*tp(:,2)-bp(:,2).*tp(:,1)];
        tdb=sum(t.*b,2);
        tdbp=sum(t.*bp,2);
        tpdb=sum(tp.*b,2);
        tpdbp=sum(tp.*bp,2);
        txtpdb=sum(txtp.*b,2);
        tpxtdbp=sum(tpxt.*bp,2);
        txbpdtp=tpxtdbp;
        tpxbdt=txtpdb;
        
        bpxtpdb=sum(bpxtp.*b,2);
        bxtdbp=sum(bxt.*bp,2);
        txbpdb=bxtdbp;
        tpxbdbp=bpxtpdb;
        txtpxt=tp-[ c c c ].*t;
        tpxtxtp=t-[ c c c ].*tp;
        txtpxbp=[ tdbp tdbp tdbp ].*tp-[tpdbp tpdbp tpdbp ].*t;
        tpxtxb=[tpdb tpdb tpdb].*t-[tdb tdb tdb].*tp;
        txbpxt=bp-[tdbp tdbp tdbp].*t;
        tpxbxtp=b-[tpdb tpdb tpdb].*tp;
        bpxtpxt=[tdbp tdbp tdbp].*tp-[c c c].*bp;
        bxtxtp=[ tpdb tpdb tpdb ].*t-[c c c].*b;
        txtpxbpxt=[tdbp tdbp tdbp].*tpxt;
        tpxtxbxtp=[tpdb tpdb tpdb].*txtp;
        txtpxbpdtp=tdbp-tpdbp.*c;
        tpxtxbdt=tpdb-tdb.*c;
        txtpxbpdb= tdbp.*tpdb-tpdbp.*tdb;
        tpxtxbdbp=txtpxbpdb;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this section calculates the coefficients for two of the forces
        temp1=tdbp.*tpdb+txtpxbpdb;
        I00a= [temp1 temp1 temp1].*tpxt;
        I00b= bxt .* [txtpxbpdtp txtpxbpdtp txtpxbpdtp];
        
        temp1=(m4pnd.*txtpdb);
        temp2=(m4pnd.*bpxtpdb);
        I_003=  [m4pd m4pd m4pd] .* I00a - [m4pnd m4pnd m4pnd] .*  I00b +  [temp1 temp1 temp1].* bpxtpxt +  [temp2 temp2 temp2].* txtpxt; 
        
        temp1=(m4pnd3.*txtpxbpdtp.*txtpdb);
        I_005=[a2m8pd a2m8pd a2m8pd] .* I00a - [a2m4pnd a2m4pnd a2m4pnd] .*  I00b - [temp1 temp1 temp1].* txtpxt;
        
        I10a=txbpxt .* [tpdb tpdb tpdb] - txtp .* [txbpdb txbpdb txbpdb];
        I10b=bxt .* [txbpdtp txbpdtp txbpdtp];
        
        temp1=(m4pn.*tdb);
        I_103= [temp1 temp1 temp1].* bpxtpxt  + m4p .* I10a -  m4pn .* I10b;
        
        temp1=m4pnd2.*(txbpdtp.*txtpdb+txtpxbpdtp.*tdb);
        I_105= a2m8p .* I10a - a2m4pn .* I10b -  [temp1 temp1 temp1].* txtpxt;
        
        I01a= txtp .* [bpxtpdb bpxtpdb bpxtpdb] - bpxtpxt.*[tpdb tpdb tpdb];
        
        temp1=(m4pn.*tpdb); 
        temp2=(m4pn.*bpxtpdb);
        I_013=   m4p .* I01a  +  [temp1 temp1 temp1] .* bpxtpxt - [temp2 temp2 temp2] .* txtp;
        
        temp1= (m4pnd2.*txtpxbpdtp.*tpdb);
        temp2=(m4pnd2.*txtpxbpdtp.*txtpdb);
        I_015= a2m8p .* I01a  -  [temp1 temp1 temp1] .* txtpxt  +  [temp2 temp2 temp2] .* txtp;
      
        temp1=(m4pnd.*txbpdtp.*tdb);
        I_205=-[temp1 temp1 temp1] .* txtpxt;
        
        temp1= (m4pnd.*txtpxbpdtp.*tpdb) ;
        I_025= [temp1 temp1 temp1].* txtp; 
        
        temp1=(m4pnd.*(txtpxbpdtp.*tdb+txbpdtp.*txtpdb));
        temp2=(m4pnd.*txbpdtp.*tpdb);
        I_115=     [temp1 temp1 temp1].* txtp  - [temp2 temp2 temp2].*txtpxt;
        
        temp1=(m4pn.*txbpdtp.*tdb);
        I_215= [temp1 temp1 temp1].* txtp;
        temp1=(m4pn.*txbpdtp.*tpdb);
        I_125=[temp1 temp1 temp1].* txtp;
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this section calculates the first two forces
        Fint_003=Fintegrals(:,2)-y(:,1).*Fintegrals(:,1);
        Fint_103=Fintegrals(:,5)-y(:,1).*Fintegrals(:,2);
        Fint_013=Fintegrals(:,4)-y(:,1).*Fintegrals(:,3);
        Fint_005=Fintegrals(:,8)-y(:,1).*Fintegrals(:,7);
        Fint_105=Fintegrals(:,11)-y(:,1).*Fintegrals(:,8);
        Fint_015=Fintegrals(:,10)-y(:,1).*Fintegrals(:,9);
        Fint_115=Fintegrals(:,13)-y(:,1).*Fintegrals(:,10);
        Fint_205=Fintegrals(:,16)-y(:,1).*Fintegrals(:,11);
        Fint_025=Fintegrals(:,14)-y(:,1).*Fintegrals(:,12);
        Fint_215=Fintegrals(:,18)-y(:,1).*Fintegrals(:,13);
        Fint_125=Fintegrals(:,15)-y(:,1).*Fintegrals(:,14);
        f4=I_003.*[Fint_003 Fint_003 Fint_003] + I_103.*[Fint_103 Fint_103 Fint_103] + I_013.*[Fint_013 Fint_013 Fint_013]; 
        f4=f4 + I_005.*[Fint_005 Fint_005 Fint_005] + I_105.*[Fint_105 Fint_105 Fint_105] + I_015.*[Fint_015 Fint_015 Fint_015];   
        f4=f4 + I_115.*[Fint_115 Fint_115 Fint_115] + I_205.*[Fint_205 Fint_205 Fint_205] + I_025.*[Fint_025 Fint_025 Fint_025]; 
        f4=f4 + I_215.*[Fint_215 Fint_215 Fint_215] + I_125.*[Fint_125 Fint_125 Fint_125];  
        f4=f4.*[oneoverL oneoverL oneoverL] ;
            
        
        Fint_003=y(:,2).*Fintegrals(:,1)-Fintegrals(:,2);
        Fint_103=y(:,2).*Fintegrals(:,2)-Fintegrals(:,5);
        Fint_013=y(:,2).*Fintegrals(:,3)-Fintegrals(:,4);
        Fint_005=y(:,2).*Fintegrals(:,7)-Fintegrals(:,8);
        Fint_105=y(:,2).*Fintegrals(:,8)-Fintegrals(:,11);
        Fint_015=y(:,2).*Fintegrals(:,9)-Fintegrals(:,10);
        Fint_115=y(:,2).*Fintegrals(:,10)-Fintegrals(:,13);
        Fint_205=y(:,2).*Fintegrals(:,11)-Fintegrals(:,16);
        Fint_025=y(:,2).*Fintegrals(:,12)-Fintegrals(:,14);
        Fint_215=y(:,2).*Fintegrals(:,13)-Fintegrals(:,18);
        Fint_125=y(:,2).*Fintegrals(:,14)-Fintegrals(:,15);
        f3=I_003.*[Fint_003 Fint_003 Fint_003] + I_103.*[Fint_103 Fint_103 Fint_103] + I_013.*[Fint_013 Fint_013 Fint_013]; 
        f3=f3 + I_005.*[Fint_005 Fint_005 Fint_005] + I_105.*[Fint_105 Fint_105 Fint_105] + I_015.*[Fint_015 Fint_015 Fint_015];   
        f3=f3 + I_115.*[Fint_115 Fint_115 Fint_115] + I_205.*[Fint_205 Fint_205 Fint_205] + I_025.*[Fint_025 Fint_025 Fint_025]; 
        f3=f3 + I_215.*[Fint_215 Fint_215 Fint_215] + I_125.*[Fint_125 Fint_125 Fint_125];  
        f3=f3.*[oneoverL oneoverL oneoverL] ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this section calculates the coefficients for the second two forces
        temp1=tpdb.* tdbp + tpxtxbdbp;
        I00a= [temp1 temp1 temp1].*txtp;
        I00b=bpxtp.* [tpxtxbdt tpxtxbdt tpxtxbdt];
        
        temp1= m4pnd .* tpxtdbp;
        temp2=m4pnd .* bxtdbp;
        I_003=  [m4pd m4pd m4pd].*I00a - [m4pnd m4pnd m4pnd].*I00b + [temp1 temp1 temp1].*bxtxtp + [temp2 temp2 temp2].*tpxtxtp;
        
        temp1=m4pnd3 .*tpxtxbdt .* tpxtdbp;
        I_005=[a2m8pd a2m8pd a2m8pd] .* I00a - [a2m4pnd a2m4pnd a2m4pnd].*  I00b -  [temp1 temp1 temp1].* tpxtxtp; 
             
        I01a= tpxt .* [tpxbdbp tpxbdbp tpxbdbp] - tpxbxtp .*[tdbp tdbp tdbp];
        I01b=-bpxtp .*[tpxbdt tpxbdt tpxbdt];
        
        temp1=m4pn .* tpdbp;
        I_013= -[temp1 temp1 temp1] .* bxtxtp + m4p .* I01a -  m4pn .* I01b;
        
        temp1=m4pnd2 .*(tpxbdt .* tpxtdbp+tpxtxbdt .* tpdbp);
        I_015= a2m8p .* I01a - a2m4pn .* I01b  +  [temp1 temp1 temp1] .* tpxtxtp;
          
        I10a=  bxtxtp .*[tdbp tdbp tdbp] - tpxt .*[bxtdbp bxtdbp bxtdbp];
        
        temp1=m4pn.*tdbp; 
        temp2=m4pn.* bxtdbp;
        I_103=   m4p .* I10a  -  [temp1 temp1 temp1].*bxtxtp + [temp2 temp2 temp2] .* tpxt;
        
        temp1=m4pnd2 .* tpxtxbdt .* tdbp;
        temp2=m4pnd2 .* tpxtxbdt .* tpxtdbp;
        I_105= a2m8p .* I10a  +  [temp1 temp1 temp1].*tpxtxtp - [temp2 temp2 temp2] .* tpxt;
        
        temp1=(m4pnd .* tpxbdt .* tpdbp);   
        I_025=-[temp1 temp1 temp1] .* tpxtxtp;
        
        temp1=(m4pnd .* tpxtxbdt .* tdbp);
        I_205= [temp1 temp1 temp1].* tpxt;
        
        temp1=m4pnd .*( tpxtxbdt .* tpdbp + tpxbdt .* tpxtdbp );
        temp2=m4pnd .* tpxbdt .* tdbp;
        I_115=[temp1 temp1 temp1] .*tpxt - [temp2 temp2 temp2] .*tpxtxtp;
        
        temp1=(m4pn .* tpxbdt .* tpdbp);
        I_125= -[temp1 temp1 temp1].* tpxt;
        
        temp1=(m4pn .* tpxbdt .* tdbp);
        I_215= -[temp1 temp1 temp1].* tpxt;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this section calculates the second two forces
        Fint_003=Fintegrals(:,3)-z(:,2).*Fintegrals(:,1);
        Fint_103=Fintegrals(:,4)-z(:,2).*Fintegrals(:,2);
        Fint_013=Fintegrals(:,6)-z(:,2).*Fintegrals(:,3);
        Fint_005=Fintegrals(:,9)-z(:,2).*Fintegrals(:,7);
        Fint_105=Fintegrals(:,10)-z(:,2).*Fintegrals(:,8);
        Fint_015=Fintegrals(:,12)-z(:,2).*Fintegrals(:,9);
        Fint_115=Fintegrals(:,14)-z(:,2).*Fintegrals(:,10);
        Fint_205=Fintegrals(:,13)-z(:,2).*Fintegrals(:,11);
        Fint_025=Fintegrals(:,17)-z(:,2).*Fintegrals(:,12);
        Fint_215=Fintegrals(:,15)-z(:,2).*Fintegrals(:,13);
        Fint_125=Fintegrals(:,19)-z(:,2).*Fintegrals(:,14);
        f1=I_003.*[Fint_003 Fint_003 Fint_003] + I_103.*[Fint_103 Fint_103 Fint_103] + I_013.*[Fint_013 Fint_013 Fint_013]; 
        f1=f1 + I_005.*[Fint_005 Fint_005 Fint_005] + I_105.*[Fint_105 Fint_105 Fint_105] + I_015.*[Fint_015 Fint_015 Fint_015];   
        f1=f1 + I_115.*[Fint_115 Fint_115 Fint_115] + I_205.*[Fint_205 Fint_205 Fint_205] + I_025.*[Fint_025 Fint_025 Fint_025]; 
        f1=f1 + I_215.*[Fint_215 Fint_215 Fint_215] + I_125.*[Fint_125 Fint_125 Fint_125];  
        f1=f1.*[oneoverLp oneoverLp oneoverLp] ;
       
        Fint_003=z(:,1).*Fintegrals(:,1)-Fintegrals(:,3);
        Fint_103=z(:,1).*Fintegrals(:,2)-Fintegrals(:,4);
        Fint_013=z(:,1).*Fintegrals(:,3)-Fintegrals(:,6);
        Fint_005=z(:,1).*Fintegrals(:,7)-Fintegrals(:,9);
        Fint_105=z(:,1).*Fintegrals(:,8)-Fintegrals(:,10);
        Fint_015=z(:,1).*Fintegrals(:,9)-Fintegrals(:,12);
        Fint_115=z(:,1).*Fintegrals(:,10)-Fintegrals(:,14);
        Fint_205=z(:,1).*Fintegrals(:,11)-Fintegrals(:,13);
        Fint_025=z(:,1).*Fintegrals(:,12)-Fintegrals(:,17);
        Fint_215=z(:,1).*Fintegrals(:,13)-Fintegrals(:,15);
        Fint_125=z(:,1).*Fintegrals(:,14)-Fintegrals(:,19);
        f2=I_003.*[Fint_003 Fint_003 Fint_003] + I_103.*[Fint_103 Fint_103 Fint_103] + I_013.*[Fint_013 Fint_013 Fint_013]; 
        f2=f2 + I_005.*[Fint_005 Fint_005 Fint_005] + I_105.*[Fint_105 Fint_105 Fint_105] + I_015.*[Fint_015 Fint_015 Fint_015];   
        f2=f2 + I_115.*[Fint_115 Fint_115 Fint_115] + I_205.*[Fint_205 Fint_205 Fint_205] + I_025.*[Fint_025 Fint_025 Fint_025]; 
        f2=f2 + I_215.*[Fint_215 Fint_215 Fint_215] + I_125.*[Fint_125 Fint_125 Fint_125];  
        f2=f2.*[oneoverLp oneoverLp oneoverLp] ;
    end
   
    if sL > 0
        % this is the parallel case the two lines are parallel use a special lower dimensional function
        [f1sp,f2sp,f3sp,f4sp]=SpecialRemoteNodeForce(x1sp,x2sp,x3sp,x4sp,bpsp,bsp,a,mu,nu,eps);
        xL=size(x1,1);
        for k=1:sL
            index=sL+1-k;
            pos=spindex(index);
            x1=[x1(1:pos-1,:); x1sp(index,:); x1(pos:xL,:)];
            x2=[x2(1:pos-1,:); x2sp(index,:); x2(pos:xL,:)];
            x3=[x3(1:pos-1,:); x3sp(index,:); x3(pos:xL,:)];
            x4=[x4(1:pos-1,:); x4sp(index,:); x4(pos:xL,:)];
            f1=[f1(1:pos-1,:); f1sp(index,:); f1(pos:xL,:)];
            f2=[f2(1:pos-1,:); f2sp(index,:); f2(pos:xL,:)];
            f3=[f3(1:pos-1,:); f3sp(index,:); f3(pos:xL,:)];
            f4=[f4(1:pos-1,:); f4sp(index,:); f4(pos:xL,:)];
            bp=[bp(1:pos-1,:); bpsp(index,:); bp(pos:xL,:)];
            b=[b(1:pos-1,:); bsp(index,:); b(pos:xL,:)];
            xL=xL+1;
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
function [f1,f2,f3,f4]=SpecialRemoteNodeForce(x1,x2,x3,x4,bp,b,a,mu,nu,ecrit);
% this calculates the forces between dislocation nodes analytically
% this is a special subroutine used for dislocation segments that are too close to parallel to be
% calculated by the regular expression for forces
%inputs: endpoints of first dislocation segment starting at x1 ending at x2 with burgers vector bp
%        endpoints of second dislocation segment starting at x3 ending at x4 with burgers vector b
%        core paramter a
%        shear modulus mu
%        poisson ration nu
%
%outputs: f1,f2,f3,f4 is the force on nodes located at x1, x2, x3, x4 respectively
     cotanthetac=sqrt((1-ecrit*(1.01))/(ecrit*1.01));
     eps=1e-16;
     Diff=x4-x3;
     oneoverL=1./sqrt(sum(Diff.*Diff,2));
     t=Diff.*[oneoverL oneoverL oneoverL];
     
     Diff=x2-x1;
     oneoverLp=1./sqrt(sum(Diff.*Diff,2));
     tp=Diff.*[oneoverLp oneoverLp oneoverLp];
     
     c=sum(t.*tp,2);
     flipL=0;
     for i=1:size(c,1)
         if c<0
             flipL=flipL+1;
             flip(flipL)=i;
             temp=x2(i,:);
             x2(i,:)=x1(i,:);
             x1(i,:)=temp;
             tp(i,:)=-tp(i,:);
             bp(i,:)=-bp(i,:);
         end
     end
     
     
     temp=sum((x2-x1).*t,2);
     x2mod=x1+[temp temp temp].*t;
     diff=(x2-x2mod);
     magdiff=sqrt(sum(diff.*diff,2));
     temp=(0.5*cotanthetac).*[magdiff magdiff magdiff].*t;
     x1mod=x1+0.5.*diff+temp;
     x2mod=x2mod+0.5.*diff-temp;
     R=(x3-x1mod);
     Rdt=sum(R.*t,2);
     nd=R-[Rdt Rdt Rdt].*t;
     d2=sum(nd.*nd,2);
     
     r4=sum(x4.*t,2);
     r3=sum(x3.*t,2);
     s2=sum(x2mod.*t,2);
     s1=sum(x1mod.*t,2);
     
     y=[ r3,  r3,  r4,  r4];
     z=[-s1, -s2, -s1, -s2];
     
     a2=a*a;
     a2_d2 = a2+d2;
     temp=1./a2_d2;
     a2d2inv=[temp temp temp temp];
     ypz=y+z;
     ymz=y-z;
     Ra    = sqrt( [a2_d2 a2_d2 a2_d2 a2_d2] + ypz.*ypz);
     Rainv=1./Ra;
     Log_Ra_ypz=log(Ra+ypz);
        
     f_003=Ra.*a2d2inv;
     f_103=-0.5.*(Log_Ra_ypz - ymz.*Ra.*a2d2inv);
     f_013=-0.5.*(Log_Ra_ypz + ymz.*Ra.*a2d2inv);
     f_113=-Log_Ra_ypz;
     f_213=z.*Log_Ra_ypz - Ra;
     f_123=y.*Log_Ra_ypz - Ra;
        
     f_005=a2d2inv.*(2.*a2d2inv.*Ra - Rainv);
     f_105= a2d2inv.*(a2d2inv.*ymz.*Ra - y.*Rainv);
     f_015=-a2d2inv.*(a2d2inv.*ymz.*Ra + z.*Rainv);
     f_115=-a2d2inv.*ypz.*Rainv;
     f_215=  Rainv - z.*f_115;
     f_125=  Rainv - y.*f_115;
     
     mf=[1 -1 -1 1]';
     Fintegrals=[f_003*mf f_103*mf f_013*mf f_113*mf f_213*mf f_123*mf f_005*mf f_105*mf f_015*mf f_115*mf f_215*mf f_125*mf];
     %            1        2        3        4        5        6        7        8        9        10       11       12 
     m4p=0.25 * mu / pi;
     m8p=0.5 * m4p;
     m4pn=m4p / ( 1 - nu );
     a2m4pn=a2*m4pn;
     a2m8p=a2 * m8p;
     
     tdb=sum(t.*b,2);
     tdbv=[tdb tdb tdb];
     tdbp=sum(t.*bp,2);
     tdbpv=[tdbp tdbp tdbp];
     nddb=sum(nd.*b,2);
     nddbv=[nddb nddb nddb];
     bxt=[b(:,2).*t(:,3)-b(:,3).*t(:,2) , b(:,3).*t(:,1)-b(:,1).*t(:,3) , b(:,1).*t(:,2)-b(:,2).*t(:,1)];
     bpxt=[bp(:,2).*t(:,3)-bp(:,3).*t(:,2) , bp(:,3).*t(:,1)-bp(:,1).*t(:,3) , bp(:,1).*t(:,2)-bp(:,2).*t(:,1)];
     ndxt=[nd(:,2).*t(:,3)-nd(:,3).*t(:,2) , nd(:,3).*t(:,1)-nd(:,1).*t(:,3) , nd(:,1).*t(:,2)-nd(:,2).*t(:,1)];
     bpxtdb=sum(bpxt.*b,2);
     bpxtdnd=sum(bpxt.*nd,2);
     bpxtdndv=[bpxtdnd bpxtdnd bpxtdnd];
     bpxtxt=tdbpv.*t - bp;
     
     I_003=m4pn.*(nddbv.*bpxtxt + [bpxtdb bpxtdb bpxtdb].*ndxt - bpxtdndv.*bxt) - m4p.*tdbv.*tdbpv.*nd; 
     I_113= (m4pn-m4p).*tdbv.*bpxtxt;
     I_005=-a2m8p.*tdbv.*tdbpv.*nd - a2m4pn.*bpxtdndv.*bxt - m4pn.*bpxtdndv.*nddbv.*ndxt;
     I_115=-a2m8p.*tdbv.*bpxtxt - m4pn.*bpxtdndv.*tdbv.*ndxt;
     
     
     Fint_003=Fintegrals(:,2)-y(:,1).*Fintegrals(:,1);
     Fint_113=Fintegrals(:,5)-y(:,1).*Fintegrals(:,4);
     Fint_005=Fintegrals(:,8)-y(:,1).*Fintegrals(:,7);
     Fint_115=Fintegrals(:,11)-y(:,1).*Fintegrals(:,10);
     
     f4=   I_003.*[Fint_003 Fint_003 Fint_003] + I_113.*[Fint_113 Fint_113 Fint_113]; 
     f4=f4+I_005.*[Fint_005 Fint_005 Fint_005] + I_115.*[Fint_115 Fint_115 Fint_115];
     f4=f4.*[oneoverL oneoverL oneoverL];
     
     Fint_003=y(:,3).*Fintegrals(:,1)-Fintegrals(:,2);
     Fint_113=y(:,3).*Fintegrals(:,4)-Fintegrals(:,5);
     Fint_005=y(:,3).*Fintegrals(:,7)-Fintegrals(:,8);
     Fint_115=y(:,3).*Fintegrals(:,10)-Fintegrals(:,11);
     
     f3=   I_003.*[Fint_003 Fint_003 Fint_003] + I_113.*[Fint_113 Fint_113 Fint_113]; 
     f3=f3+I_005.*[Fint_005 Fint_005 Fint_005] + I_115.*[Fint_115 Fint_115 Fint_115];
     f3=f3.*[oneoverL oneoverL oneoverL];
     
     corsize=0;
     for i=1:size(c,1)
         if (diff(i,:)*diff(i,:)')>(eps*(x2mod(i,:)*x2mod(i,:)'+x1mod(i,:)*x1mod(i,:)'))
             corsize=corsize+1;
             corindex(corsize)=i;
             x1mod2(corsize,:)=x1mod(i,:);
             x12(corsize,:)=x1(i,:);
             x2mod2(corsize,:)=x2mod(i,:);
             x22(corsize,:)=x2(i,:);
             x32(corsize,:)=x3(i,:);
             x42(corsize,:)=x4(i,:);
             bp2(corsize,:)=bp(i,:);
             b2(corsize,:)=b(i,:);
         end
     end
     if corsize>0
         [whocares1,whocares2,f3cor2,f4cor2]=RemoteNodeForce(x12,x1mod2,x32,x42,bp2,b2,a,mu,nu);
         [whocares1,whocares2,f3cor3,f4cor3]=RemoteNodeForce(x2mod2,x22,x32,x42,bp2,b2,a,mu,nu);
         f3cor=zeros(size(c,1),3);
         f4cor=zeros(size(c,1),3);
         for i=1:corsize
             index=corindex(i);
             f3cor(index,:)=f3cor2(i,:)+f3cor3(i,:);
             f4cor(index,:)=f4cor2(i,:)+f4cor3(i,:);
         end
         f3=f3+f3cor;
         f4=f4+f4cor;
     end
     
     
     
     temp=sum((x4-x3).*tp,2);
     x4mod=x3+[temp temp temp].*tp;
     diff=x4-x4mod;
     magdiff=sqrt(sum(diff.*diff,2));
     temp=(0.5*cotanthetac).*[magdiff magdiff magdiff].*tp;
     x3mod=x3+0.5.*diff+temp;
     x4mod=x4mod+0.5.*diff-temp;
     R=(x3mod-x1);
     Rdtp=sum(R.*tp,2);
     nd=R-[Rdtp Rdtp Rdtp].*tp;
     d2=sum(nd.*nd,2);
     r4=sum(x4mod.*tp,2);
     r3=sum(x3mod.*tp,2);
     s2=sum(x2.*tp,2);
     s1=sum(x1.*tp,2);
     
     y=[ r3,  r3,  r4,  r4];
     z=[-s1, -s2, -s1, -s2];
     
     a2=a*a;
     a2_d2 = a2+d2;
     temp=1./a2_d2;
     a2d2inv=[temp temp temp temp];
     ypz=y+z;
     ymz=y-z;
     Ra    = sqrt( [a2_d2 a2_d2 a2_d2 a2_d2] + ypz.*ypz);
     Rainv=1./Ra;
     Log_Ra_ypz=log(Ra+ypz);
        
     f_003=Ra.*a2d2inv;
     f_103=-0.5.*(Log_Ra_ypz - ymz.*Ra.*a2d2inv);
     f_013=-0.5.*(Log_Ra_ypz + ymz.*Ra.*a2d2inv);
     f_113=-Log_Ra_ypz;
     f_213=z.*Log_Ra_ypz - Ra;
     f_123=y.*Log_Ra_ypz - Ra;
        
     f_005=a2d2inv.*(2.*a2d2inv.*Ra - Rainv);
     f_105= a2d2inv.*(a2d2inv.*ymz.*Ra - y.*Rainv);
     f_015=-a2d2inv.*(a2d2inv.*ymz.*Ra + z.*Rainv);
     f_115=-a2d2inv.*ypz.*Rainv;
     f_215=  Rainv - z.*f_115;
     f_125=  Rainv - y.*f_115;
     
     mf=[1 -1 -1 1]';
     Fintegrals=[f_003*mf f_103*mf f_013*mf f_113*mf f_213*mf f_123*mf f_005*mf f_105*mf f_015*mf f_115*mf f_215*mf f_125*mf];
     %            1        2        3        4        5        6        7        8        9        10       11       12 
     
     tpdb=sum(tp.*b,2);
     tpdbv=[tpdb tpdb tpdb];
     tpdbp=sum(tp.*bp,2);
     tpdbpv=[tpdbp tpdbp tpdbp];
     nddbp=sum(nd.*bp,2);
     nddbpv=[nddbp nddbp nddbp];
     bxtp=[b(:,2).*tp(:,3)-b(:,3).*tp(:,2) , b(:,3).*tp(:,1)-b(:,1).*tp(:,3) , b(:,1).*tp(:,2)-b(:,2).*tp(:,1)];
     bpxtp=[bp(:,2).*tp(:,3)-bp(:,3).*tp(:,2) , bp(:,3).*tp(:,1)-bp(:,1).*tp(:,3) , bp(:,1).*tp(:,2)-bp(:,2).*tp(:,1)];
     ndxtp=[nd(:,2).*tp(:,3)-nd(:,3).*tp(:,2) , nd(:,3).*tp(:,1)-nd(:,1).*tp(:,3) , nd(:,1).*tp(:,2)-nd(:,2).*tp(:,1)];
     bxtpdbp=sum(bxtp.*bp,2);
     bxtpdnd=sum(bxtp.*nd,2);
     bxtpdndv=[bxtpdnd bxtpdnd bxtpdnd];
     bxtpxtp=tpdbv.*tp - b;
     
     I_003=m4pn.*(nddbpv.*bxtpxtp + [bxtpdbp bxtpdbp bxtpdbp].*ndxtp - bxtpdndv.*bpxtp) - m4p.*tpdbpv.*tpdbv.*nd; 
     I_113= (m4pn-m4p).*tpdbpv.*bxtpxtp;
     I_005=-a2m8p.*tpdbpv.*tpdbv.*nd - a2m4pn.*bxtpdndv.*bpxtp - m4pn.*bxtpdndv.*nddbpv.*ndxtp;
     I_115=-a2m8p.*tpdbpv.*bxtpxtp - m4pn.*bxtpdndv.*tpdbpv.*ndxtp;
     
     
     Fint_003=Fintegrals(:,3)-z(:,1).*Fintegrals(:,1);
     Fint_113=Fintegrals(:,6)-z(:,1).*Fintegrals(:,4);
     Fint_005=Fintegrals(:,9)-z(:,1).*Fintegrals(:,7);
     Fint_115=Fintegrals(:,12)-z(:,1).*Fintegrals(:,10);
     
     f2=   I_003.*[Fint_003 Fint_003 Fint_003] + I_113.*[Fint_113 Fint_113 Fint_113]; 
     f2=f2+I_005.*[Fint_005 Fint_005 Fint_005] + I_115.*[Fint_115 Fint_115 Fint_115];
     f2=f2.*[oneoverLp oneoverLp oneoverLp];
     
     Fint_003=z(:,2).*Fintegrals(:,1)-Fintegrals(:,3);
     Fint_113=z(:,2).*Fintegrals(:,4)-Fintegrals(:,6);
     Fint_005=z(:,2).*Fintegrals(:,7)-Fintegrals(:,9);
     Fint_115=z(:,2).*Fintegrals(:,10)-Fintegrals(:,12);
     
     f1=   I_003.*[Fint_003 Fint_003 Fint_003] + I_113.*[Fint_113 Fint_113 Fint_113]; 
     f1=f1+I_005.*[Fint_005 Fint_005 Fint_005] + I_115.*[Fint_115 Fint_115 Fint_115];
     f1=f1.*[oneoverLp oneoverLp oneoverLp];
     
     corsize=0;
     corindex=[];
     x12=[];
     x22=[];
     x42=[];
     x4mod2=[];
     x3mod2=[];
     x32=[];
     bp2=[];
     b2=[];
     for i=1:size(c,1)
         if (diff(i,:)*diff(i,:)')>(eps*(x4mod(i,:)*x4mod(i,:)'+x3mod(i,:)*x3mod(i,:)'))
             corsize=corsize+1;
             corindex(corsize)=i;
             x12(corsize,:)=x1(i,:);
             x22(corsize,:)=x2(i,:);
             x32(corsize,:)=x3(i,:);
             x3mod2(corsize,:)=x3mod(i,:);
             x42(corsize,:)=x4(i,:);
             x4mod2(corsize,:)=x4mod(i,:);
             bp2(corsize,:)=bp(i,:);
             b2(corsize,:)=b(i,:);
         end
     end
     if corsize>0
         [whocares1,whocares2,f1cor2,f2cor2]=RemoteNodeForce(x32,x3mod2,x12,x22,b2,bp2,a,mu,nu);
         [whocares1,whocares2,f1cor3,f2cor3]=RemoteNodeForce(x4mod2,x42,x12,x22,b2,bp2,a,mu,nu);
         f1cor=zeros(size(c,1),3);
         f2cor=zeros(size(c,1),3);
         for i=1:corsize
             index=corindex(i);
             f1cor(index,:)=f1cor2(i,:)+f1cor3(i,:);
             f2cor(index,:)=f2cor2(i,:)+f2cor3(i,:);
         end
         f1=f1+f1cor;
         f2=f2+f2cor;
     end
     
     for i=1:flipL
         index=flip(i);
         temp=x2(index,:);
         x2(index,:)=x1(index,:);
         x1(index,:)=temp;
         bp(index,:)=-bp(index,:);
         temp=f2(index,:);
         f2(index,:)=f1(index,:);
         f1(index,:)=temp;
     end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     