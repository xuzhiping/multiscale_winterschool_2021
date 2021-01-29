clear all;
close all;
tic;
% create result directory
if ~exist('./Result','dir')
    mkdir ./Result
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=256;  % 
nx=n;
ny=n;
x=linspace(1,n,n);  % 
y=linspace(1,n,n);  % 
xs=x(2)-x(1);     % 
ys=y(2)-y(1);     % 
fx=2*pi/xs;    %
fy=2*pi/ys;    %
kx=fx*(1:n)/n;  %
ky=fy*(1:n)/n;  %
vect=ones(n);
k2x=zeros(n,n);
k2y=zeros(n,n);
   for i=1:1:n; 
   k2x(i,:)=kx(i)*vect(i,:);
   end
for i=1:1:n;
   k2y(:,i)=ky(i)*vect(:,i);
end
 e=1e-30;
 kkx=(-2*(cos(k2x)-1))+e/2;  % 
 kky=(-2*(cos(k2y)-1))+e/2;  % 
 kk=kkx+kky;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some related parameters
timetot = 1000;
pnum_s=4;
p_every=20;
step = 0.05;
mobile = 1;
raodong = 0.001;
beta = 1;
alfa1 = 1;
alfa2 = 1;
alfa3 = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 10; % 
eta=zeros(nx,ny,p);
chem_eta=zeros(nx,ny,p);
f_chem_eta=zeros(nx,ny,p);
f_eta=zeros(nx,ny,p);
mid_feta=zeros(nx,ny,p);
nframe=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:1:p;
%     eta(:,:,i) = zeros(nx,ny);% + 0.005*(0.5-rand(nx,ny));
% end

 for m=1:1:p; 
ax=round(rand(1,1)*nx);
ay=round(rand(1,1)*ny);  
   for i = 1:1:nx;
       for j=1:1:ny;

        if round(sqrt((i-ax)^2+(j-ay)^2)) <= 15               
                     eta(i,j,m)=1.0;
        end
       end  
   end
 end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 s_eta2 = 0; 
 for m=1:1:p;
   s_eta2 = s_eta2 + eta(:,:,m);
 end 
% h=figure;
% pcolor(s_eta2);    
% shading interp;
% colormap(jet);
% colorbar

flag1=s_eta2>=1.1; 
flag2(:,:)=flag1;
 for m=1:1:p;
   eta(:,:,m) = eta(:,:,m).*(1-flag1);
 end 
 s_eta2 = 0; 
 for m=1:1:p;
   s_eta2 = s_eta2 + eta(:,:,m);
 end 
% h=figure;
% pcolor(s_eta2);    
% shading interp;
% colormap(jet);
% colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for j=1:1:timetot;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
     s_eta2 = 0; 
     for m=1:1:p;
       s_eta2 = s_eta2 + eta(:,:,m).^2;
     end    
    
     for g = 1:1:p;
         chem_eta(:,:,g) = -alfa1*eta(:,:,g)+alfa2*eta(:,:,g).^3 + 2*alfa3*eta(:,:,g).*(s_eta2-eta(:,:,g).^2);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
         f_chem_eta(:,:,g)=fft2(chem_eta(:,:,g));
         f_eta(:,:,g)=fft2(eta(:,:,g));
         mid_feta(:,:,g)=(f_eta(:,:,g)-step*mobile.*f_chem_eta(:,:,g))./(1+step*mobile*beta*kk); 
         eta(:,:,g)=real(ifft2(mid_feta(:,:,g)));%+raodong*(0.5-rand(nx,ny));
     end

   if mod(j,timetot/pnum_s) == 0
        num = j/(timetot/pnum_s)
%         h=figure;
%         pcolor(s_eta2);    
%         shading interp;
%         colormap(jet);
%         colorbar
        
        a=s_eta2;
        b = mat2gray(a,[min(min(a)) max(max(a))]); % Convert a matrix to a grayscale intensity image
        dd = gray2ind(b); % Convert an intensity image to an indexed image       
        name = ['./Result/profile' num2str(num) '.jpg'];% 
        imwrite(dd,jet,name,'jpg') % Write two-dimensional matrix to graphics file 
        
   end
     if mod(j,p_every) ==0;
            nframe=nframe+1;
            a = s_eta2;
            b = mat2gray(a,[min(min(a)) max(max(a))]); % Convert a matrix to a grayscale intensity image
             dd = gray2ind(b); % Convert an intensity image to an indexed image  
             imshow(dd,jet,'border','tight');
             frame=getframe(gcf);
             im{nframe}=frame2im(frame);
     end    
end
for idx = 1:nframe
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,['./Result/animate.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
    else
        imwrite(A,map,['./Result/animate.gif'],'gif','WriteMode','append','DelayTime',0.5);
    end
end

toc