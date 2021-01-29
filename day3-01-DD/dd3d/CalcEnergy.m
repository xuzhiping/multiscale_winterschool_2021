function [W]=CalcEnergy(rn,links,mu,nu,a,Ec);
% this function calculates the internal energy stored in the dislocation network
% the inputs are rn, links, mu, nu, a, and Ec
% the output is energy W
% function [W]=CalcEnergy(rn,links,mu,nu,a,Ec)

    L=size(links,1);
    Wself=SegSelfEnergy(mu,nu,Ec,a,rn(links(:,1),1:3),rn(links(:,2),1:3),links(:,3:5));
    index1=[];index2=[];index3=[];index4=[];bindex=[];bpindex=[];
    for i=1:L-1
        index1=[index1;ones(L-i,1)*links(i,1)];  index2=[index2;ones(L-i,1)*links(i,2)];
        index3=[index3;links((i+1):L,1)];        index4=[index4;links((i+1):L,2)];
        bpindex=[bpindex;ones(L-i,1)*i];         bindex=[bindex;(linspace(i+1,L,L-i))'];    
    end
    Winteract=SegSegInteractionEnergy(rn(index1,1:3),rn(index2,1:3),rn(index3,1:3),...
        rn(index4,1:3),links(bpindex,3:5),links(bindex,3:5),a,mu,nu);
    W=sum(Wself)+sum(Winteract);



    
function [W]=SegSelfEnergy(mu,nu,Ec,a,x1,x2,b);
% this function is called by CalcEnergy and it calculates the self energy of a dislocation segment
% inputs are mu, nu, Ec, a, x1, x2, and b
% x1, x2, and b are three matrices made up of x1, x2, b row vector triplets
% the output is a vector W that has a length of the number of triplets
% function [W]=SegSelfEnergy(mu,nu,Ec,a,x1,x2,b)

    Diff=x2-x1;

    L=sqrt(sum(Diff.*Diff,2));
    Linv=1./L;
    t=Diff.*[Linv Linv Linv];

    bs=sum(b.*t,2);
    bs2=bs.*bs;
    be2=sum(b.*b,2)-bs2;

    La=sqrt(L.*L+a*a);
    common=2.*L.*log((La+L)./a);
    %Self Elastic Energy
    W=(0.125*mu/pi) .*(  be2.*common./(1-nu) + bs2.*(common -((3-nu)/(1-nu)).*(La - a) ));
    %Self Core energy
    W=W+(bs2 + be2./(1-nu)).*Ec.*L;
    
    

    
function [W12]=SegSegInteractionEnergy(x1,x2,x3,x4,bp,b,a,mu,nu);
% this subroutine calculates the interaction energy between two dislocation line segments
% segment1 starts at x1 and ends at x2 with burgers vector bp 
% segment2 starts at x3 and ends at x4 with burgers vector b
% a is the core width
% mu is the shear modulus
% nu is poisson's ratio
% the function is used with matrices of row vectors
% the output is given in W12
%[W12]=SegSegInteractionEnergy(x1,x2,x3,x4,bp,b,a,mu,nu);

    W12=zeros(size(x1,1),1);
    eps=1e-6;
    
    Diff=x4-x3;
    oneoverL=1./sqrt(sum(Diff.*Diff,2));
    t=Diff.*[oneoverL oneoverL oneoverL];
    Diff=x2-x1;
    oneoverLp=1./sqrt(sum(Diff.*Diff,2));
    tp=Diff.*[oneoverLp oneoverLp oneoverLp];

    c=sum(t.*tp,2);
    onemc2=1-c.*c;  
    cL=size(c,1);
    k=1;
    spindex=[];
    for i=1:cL
        index=cL+1-i;
        if onemc2(index)<eps
            spindex(k)=index;
            
            x1sp(k,:)=x1(index,:); x2sp(k,:)=x2(index,:); x3sp(k,:)=x3(index,:); x4sp(k,:)=x4(index,:);
            bpsp(k,:)=bp(index,:); bsp(k,:)=b(index,:);
            
            x1(index,:)=[]; x2(index,:)=[]; x3(index,:)=[]; x4(index,:)=[];
            bp(index,:)=[]; b(index,:)=[];
            
            t(index,:)=[]; tp(index,:)=[];
            c(index,:)=[]; onemc2(index,:)=[];
            
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
        a2_d2v=[a2_d2 a2_d2 a2_d2 a2_d2];
        cv=[c c c c];
        onemc2inv2=[onemc2inv onemc2inv onemc2inv onemc2inv];
        
        Ra    = sqrt( a2_d2v + yin.*yin + zin.*zin + 2.*yin.*zin.*cv );        
        Rainv = 1./Ra;      
         
         log_Ra_Rdot_tp =      log(Ra+zin+yin.*cv);
        ylog_Ra_Rdot_tp = yin.*log_Ra_Rdot_tp;
        
         log_Ra_Rdot_t  =      log(Ra+yin+zin.*cv);
        zlog_Ra_Rdot_t  = zin.*log_Ra_Rdot_t;
       
        denom=1./sqrt(onemc2.*a2_d2);  
        cdenom=(1+c).*denom;
        
        W_003=-2.*[denom denom denom denom].*atan((Ra+yin+zin).*[cdenom cdenom cdenom cdenom]);
        
        adW_003=a2_d2v.*W_003;
        commonW223=(cv.*Ra - adW_003 ).*onemc2inv2;
        W_001=zlog_Ra_Rdot_t +  ylog_Ra_Rdot_tp - adW_003;
        W_103=( cv.*log_Ra_Rdot_t  - log_Ra_Rdot_tp ).*onemc2inv2;
        W_013=( cv.*log_Ra_Rdot_tp - log_Ra_Rdot_t  ).*onemc2inv2;
        W_113=( cv.*adW_003 - Ra ).*onemc2inv2;
        W_203= zlog_Ra_Rdot_t  + commonW223;
        W_023= ylog_Ra_Rdot_tp + commonW223;
        
        mf=[1 -1 -1 1]';
        Wintegrals=[W_001*mf W_003*mf W_103*mf W_013*mf W_113*mf W_203*mf W_023*mf ];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this section calculates the dot products and cross prodcucts for the coefficients
        
        omn=1-nu;
        a2omnd2=0.5*a2*omn;
     
        tdb =sum(t.*b,2);
        tdbp=sum(t.*bp,2);
        tpdb =sum(tp.*b,2);
        tpdbp=sum(tp.*bp,2);
        bdbp =sum(b.*bp,2);
        txtpdb =sum(txtp.*b,2);
        txtpdbp=sum(txtp.*bp,2);
        
        common1=tdb.*tpdbp;
        common2=tdbp.*tpdb;
        
        I_001=omn.*common1+(2*nu).*common2-bdbp.*c;
        I_003=a2omnd2.*common1+(a2*nu).*common2-(a2.*bdbp+d.*d.*txtpdb.*txtpdbp).*c;
        I_103=-d.*c.*(tdb  .* txtpdbp + tdbp .* txtpdb);
        I_013=-d.*c.*(tpdb .* txtpdbp + tpdbp.* txtpdb);
        I_113=-c.*(common2+common1);
        I_203=-tdb.*tdbp.*c;
        I_023=-tpdb.*tpdbp.*c;
        
        W12=I_001.*Wintegrals(:,1)+I_003.*Wintegrals(:,2)+I_103.*Wintegrals(:,3)+I_013.*Wintegrals(:,4);
        W12=W12+I_113.*Wintegrals(:,5)+I_203.*Wintegrals(:,6)+I_023.*Wintegrals(:,7);
        W12=W12.*(0.25*mu/(pi*(1-nu)));
    end
   
    if sL > 0
        % this is the parallel case the two lines are parallel use a special lower dimensional function
        [W12sp]=SpecialSegSegInteractionEnergy(x1sp,x2sp,x3sp,x4sp,bpsp,bsp,a,mu,nu,eps);
        xL=size(x1,1);
        for k=1:sL
            index=sL+1-k;
            pos=spindex(index);
            x1=[x1(1:pos-1,:); x1sp(index,:); x1(pos:xL,:)];
            x2=[x2(1:pos-1,:); x2sp(index,:); x2(pos:xL,:)];
            x3=[x3(1:pos-1,:); x3sp(index,:); x3(pos:xL,:)];
            x4=[x4(1:pos-1,:); x4sp(index,:); x4(pos:xL,:)];
            W12=[W12(1:pos-1); W12sp(index); W12(pos:xL)];
            bp=[bp(1:pos-1,:); bpsp(index,:); bp(pos:xL,:)];
            b=[b(1:pos-1,:); bsp(index,:); b(pos:xL,:)];
            xL=xL+1;
        end
    end
    
    
    
    
function [W12]=SpecialSegSegInteractionEnergy(x1,x2,x3,x4,bp,b,a,mu,nu,ecrit);
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
     x1mod=x1;
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
     Ra    = sqrt( [a2_d2 a2_d2 a2_d2 a2_d2] + ypz.*ypz);
     Log_Ra_ypz=log(Ra+ypz);
        
     f_001=ypz.*Log_Ra_ypz-Ra;
     f_003=Ra.*a2d2inv;
     f_113=-Log_Ra_ypz;
     f_223=ypz.*Log_Ra_ypz-2.*Ra;
     
     mf=[1 -1 -1 1]';
     Wintegrals=[f_001*mf f_003*mf f_113*mf f_223*mf];
     %            1        2        3        4       
     opn=1+nu;
     a2opnd2=0.5*a2*opn;
     
     tdb=sum(t.*b,2);
     tdbp=sum(t.*bp,2);
     bdbp=sum(b.*bp,2);
     nddb=sum(nd.*b,2);
     nddbp=sum(nd.*bp,2);
     
     I_001=opn.*tdb.*tdbp-bdbp;
     I_003=a2opnd2.*tdb.*tdbp - a2*bdbp - nddb.*nddbp;
     I_113=-tdb.*nddbp-tdbp.*nddb;
     I_223=-tdb.*tdbp;
     W12=(0.25*mu/(pi*(1-nu))).*(I_001.*Wintegrals(:,1)+I_003.*Wintegrals(:,2)+I_113.*Wintegrals(:,3)+I_223.*Wintegrals(:,4));
     
     corsize=0;
     for i=1:size(c,1)
         if (diff(i,:)*diff(i,:)')>(eps*(x2mod(i,:)-x1mod(i,:))*(x2mod(i,:)-x1mod(i,:))')
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
        W12corr=SegSegInteractionEnergy(x12,x1mod2,x32,x42,bp2,b2,a,mu,nu);
        W12corr=W12corr+SegSegInteractionEnergy(x2mod2,x22,x32,x42,bp2,b2,a,mu,nu);
        for i=1:corsize
            index=corindex(i);
            W12(index)=W12(index)+W12corr(i);
        end
     end
     
     for i=1:flipL
         index=flip(i);
         temp=x2(index,:);
         x2(index,:)=x1(index,:);
         x1(index,:)=temp;
         bp(index,:)=-bp(index,:);
     end     