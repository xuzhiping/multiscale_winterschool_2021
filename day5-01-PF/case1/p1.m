clear all;
close all;
tic;
 %设置控制方程计算参数
wa=2.3;
beta=1;
unitt=0.02;  %
timetot=10000; 
pnum_s=4;
p_every=500;
% 创建结果文件夹
if ~exist('./Result','dir')
    mkdir ./Result
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%离散实空间和傅里叶空间
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=256;   
h=n;   
x=linspace(1,h,n);   
y=linspace(1,h,n);   
xs=x(2)-x(1);      
ys=y(2)-y(1);      
fx=2*pi/xs;    
fy=2*pi/ys;    
kx=fx*(1:n)/n;  
ky=fy*(1:n)/n;  
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
 kkx=(-2*(cos(k2x)-1))+e/2;  
 kky=(-2*(cos(k2y)-1))+e/2;  
 kk=kkx+kky;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %设置初值
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
turbulence=0.01; % 
c0=0.5; % 
c=c0+turbulence*(0.5-rand(n));
c1=c;
nframe=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ics=1:1:timetot;

   p= wa*(1-2*c)+log(complex(c)./(1-complex(c)));
   fc=fft2(c);  % 
   fp=fft2(p);  
   midc = (fc + unitt*(-kk.*fp-beta*kk.^2.*fc));
   c=ifft2(midc); 
   c=real(c);
%    c=c-mean(mean(c-c1)); % 
   
    if ics<=100
        c=c+turbulence*(randn(n));
    end
   
     if mod(ics,timetot/pnum_s) ==0;
            j = ics/(timetot/pnum_s)
            a = gather(c);
            b = mat2gray(a,[min(min(a)) max(max(a))]); % Convert a matrix to a grayscale intensity image
             dd = gray2ind(b); % Convert an intensity image to an indexed image       
             name = ['./Result/profile' num2str(j) '.jpg'];% 
             imwrite(dd,jet,name,'jpg') % Write two-dimensional matrix to graphics file  
     end 
     
     if mod(ics,p_every) ==0;
            nframe=nframe+1;
            a = gather(c);
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







