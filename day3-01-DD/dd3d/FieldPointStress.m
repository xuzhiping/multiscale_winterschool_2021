function [Stress]=FieldPointStress(x,x1,x2,b,a,mu,nu);
% inputs:x(dim,1:3) a list of positions where you want the stress calculated
%        x1(nseg,1:3) a list of starting positions of the dislocations segments
%        x2(nseg,1:3) a list of ending positions of the dislocation segments
%        b(nseg,1:3) a list of burgers vectors associated with the segments going from x1 to x2
%        a=core width
%        mu= shear modulus
%        nu=poisson's ratio
%
% outputs:Stress(dim,1:6) stresses at the field points requested it will calculate the field point stress as a sum of all the segments
% the components of the stress are organized in the following fashion
% Stress(i,:)=[s_11 s_22 s_33 s_12 s_23 s_13]
sx=size(x,1);
if sx>1
     for i=1:sx
         [Stress(i,:)]=FieldPointStress(x(i,:),x1,x2,b,a,mu,nu);
     end
else
     sx1=size(x1,1);
     x=ones(sx1,1)*x;   

     Diff=x2-x1;
     oneoverL=1./sqrt(sum(Diff.*Diff,2));
     t=Diff.*[oneoverL oneoverL oneoverL];
     
     R=(x-x1);
     Rdt=sum(R.*t,2);
     nd=R-[Rdt Rdt Rdt].*t;
     d2=sum(nd.*nd,2);
     
     s(:,1)=-sum(R.*t,2);
     s(:,2)=-sum((x-x2).*t,2);
    
     a2=a*a;
     a2_d2 = a2+d2;
     temp=1./a2_d2;
     a2d2inv=[temp temp];
     
     Ra    = sqrt( [a2_d2 a2_d2] + s.*s);
     Rainv=1./Ra;
     Ra3inv=Rainv.*Rainv.*Rainv;
     sRa3inv=s.*Ra3inv;
     
     s_03=s.*Rainv.*a2d2inv;
     s_13=-Rainv;
     s_05=(2.*s_03+sRa3inv).*a2d2inv;
     s_15=-Ra3inv;
     s_25=s_03-sRa3inv;
     
     mf=[-1 1]';
     s_03=s_03*mf;
     s_13=s_13*mf;
     s_05=s_05*mf;
     s_15=s_15*mf;
     s_25=s_25*mf;
     
     m4p=0.25 * mu / pi;
     m8p=0.5 * m4p;
     m4pn=m4p / ( 1 - nu );
     mn4pn=m4pn*nu;
     a2m4pn=a2 * m4pn;
     a2m8p=a2 * m8p;
     
     onev=ones(size(x1,1),1);
     zerov=zeros(size(x1,1),1);
     eyesix=[onev onev onev zerov zerov zerov];
     
     txb=[ t(:,2).*b(:,3)-t(:,3).*b(:,2) ,   t(:,3).*b(:,1)-t(:,1).*b(:,3) ,   t(:,1).*b(:,2)-t(:,2).*b(:,1)];
     dxb=[nd(:,2).*b(:,3)-nd(:,3).*b(:,2) , nd(:,3).*b(:,1)-nd(:,1).*b(:,3) , nd(:,1).*b(:,2)-nd(:,2).*b(:,1)];
     dxbdt=sum(dxb.*t,2);
     dxbdtv=[dxbdt dxbdt dxbdt dxbdt dxbdt dxbdt];
     
     dmd=[     nd(:,1).* nd(:,1)    nd(:,2).* nd(:,2)     nd(:,3).* nd(:,3) nd(:,1).* nd(:,2)                     nd(:,2).* nd(:,3)                   nd(:,3).* nd(:,1)];
     tmt=[      t(:,1).*  t(:,1)     t(:,2).*  t(:,2)      t(:,3).*  t(:,3)  t(:,1).*  t(:,2)                      t(:,2).*  t(:,3)                    t(:,3).*  t(:,1)];
     tmd=[   2.*t(:,1).* nd(:,1)  2.*t(:,2).* nd(:,2)   2.*t(:,3).* nd(:,3)  t(:,1).* nd(:,2)+  nd(:,1).* t(:,2)    t(:,2).* nd(:,3)+  nd(:,2).* t(:,3)  t(:,3).* nd(:,1)+ nd(:,3).* t(:,1)];
     tmtxb=[ 2.*t(:,1).*txb(:,1)  2.*t(:,2).*txb(:,2)   2.*t(:,3).*txb(:,3)  t(:,1).*txb(:,2)+ txb(:,1).* t(:,2)    t(:,2).*txb(:,3)+ txb(:,2).* t(:,3)  t(:,3).*txb(:,1)+txb(:,3).* t(:,1)];
     dmtxb=[2.*nd(:,1).*txb(:,1) 2.*nd(:,2).*txb(:,2)  2.*nd(:,3).*txb(:,3) nd(:,1).*txb(:,2)+ txb(:,1).*nd(:,2)   nd(:,2).*txb(:,3)+ txb(:,2).*nd(:,3) nd(:,3).*txb(:,1)+txb(:,3).*nd(:,1)];
     tmdxb=[ 2.*t(:,1).*dxb(:,1)  2.*t(:,2).*dxb(:,2)   2.*t(:,3).*dxb(:,3)  t(:,1).*dxb(:,2)+ dxb(:,1).* t(:,2)    t(:,2).*dxb(:,3)+ dxb(:,2).* t(:,3)  t(:,3).*dxb(:,1)+dxb(:,3).* t(:,1)];
     
     common=m4pn.*dxbdtv;
     
     I_03=m4pn.*(dxbdtv.*eyesix+dmtxb)-m4p.*tmdxb;
     I_13=-mn4pn.*tmtxb;
     I_05=common.*(a2.*eyesix+dmd)-a2m8p.*tmdxb;
     I_15=a2m8p.*tmtxb-common.*tmd;
     I_25=common.*tmt;
     
     Stress=       I_03.*[s_03 s_03 s_03 s_03 s_03 s_03]+I_13.*[s_13 s_13 s_13 s_13 s_13 s_13];
     Stress=Stress+I_05.*[s_05 s_05 s_05 s_05 s_05 s_05]+I_15.*[s_15 s_15 s_15 s_15 s_15 s_15]+I_25.*[s_25 s_25 s_25 s_25 s_25 s_25];
     Stress=sum(Stress,1);
 end