%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
tic;
% create result directory
if ~exist('./Result','dir')
    mkdir ./Result
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=1024; 
nx = n; % the total number of grids along direction x
ny = n; % the total number of grids along direction y
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

 kkx=(-2*(cos(k2x)-1));  
 kky=(-2*(cos(k2y)-1)); 
 kk=kkx+kky;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 timetot = 1e3;
 pnum_s=8;
 step = 0.5;
 epsilon = 0.15;
 lro0 = -0.23;
 lro = (lro0*ones(nx,ny));
 change = 1e-2;
 p_every=20;
 nframe=0;

 for i=1:1:timetot;
     f_lro=fft2(lro);
     mid_flro = (f_lro-step*kk.*fft2(lro.^3))./(1-step*(-kk.^5+14/3*kk.^4-73/9*kk.^3+56/9*kk.^2-(16/9-epsilon)*kk));
     lro=real(ifft2(mid_flro))+change*(0.5-rand(nx,ny));
     
     if mod(i,timetot/pnum_s) ==0;
                        
            j = i/(timetot/pnum_s)
%             h=figure;
%             pcolor(lro);    
%             shading interp;
%             colormap(jet);
%             colorbar        

             a=lro;
             b = mat2gray(a,[min(min(a)) max(max(a))]); % Convert a matrix to a grayscale intensity image
             dd = gray2ind(b); % Convert an intensity image to an indexed image       
             name = ['./Result/profile' num2str(j) '.jpg'];% 
             imwrite(dd,jet,name,'jpg') % Write two-dimensional matrix to graphics file  
     end     
      if mod(i,p_every) ==0;
            nframe=nframe+1;
            a = lro;
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
     