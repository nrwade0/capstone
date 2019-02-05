%/
% Step 2: rotate data 
% input: data from pa_fit and pix size
% receives and subsequently rotates velocity and intensity maps of NGC 4736
% Takes vel and int maps to rotate so the major kinematic
% axis is along the positive x-axis.  This rotates the data by the pa.
% output: two .fits maps of rotated intensity and velocity maps.
% also outputs goodbye=1 to ensure completion 
%/

clear
clc
format long g

% load step 1 data
% on home pc
load('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\matlab_data\pa_fit.mat')
% on mlc mac w/ usb
% load('/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/matlab_data/pa_fit.mat');

% correct velocity data by flipping up/down
% onto home pc
v=flipud(fitsread('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\fits_data\n4321_v.fits'));
% on mlc mac w/ usb
% v=flipud(fitsread('/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/fits_data/n4321_v.fits'));

% rid any of any values below 0 (=nan)
v(v==0)=nan;
% correct for inclination angle and systemic velocity
vy=(v-vsys)./sind(ia);

% intensity data input (corrected for flip)
% on home pc
I=flipud(fitsread('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\fits_data\n4321_i.fits'));
% on mlc mac w/ usb
% I=flipud(fitsread('/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/fits_data/n4321_i.fits'));

% pix size from DS9
pix_size=0.000223611*3600;

% check data to make sure it makes sense w/ DS9
%
% %velocity data corrected
% figure(1)
% clf
% imagesc(vy)
% colorbar
% colormap('jet')
%
% % intensity data corrected (log scale, abs value)
% figure(2)
% clf
% imagesc(log(abs(I)))
% colorbar
% colormap('jet')

[row,col]=size(I);

x=zeros(row,col); 
y=zeros(row,col);

for i=1:col
    x(:,i)=i;  
end
for i=1:row
    y(i,:)=row-i+1; 
end

% create X, Y coords in galaxy plane
X=(-(x-xo).*sind(pa)+(y-yo).*cosd(pa)).*pix_size;
Y=((-(x-xo).*cosd(pa)-(y-yo).*sind(pa)))./cosd(ia).*pix_size;
R=sqrt(X.^2+Y.^2);

% % gradient from top left to bottom right proves what?
% figure(4)
% clf
% imagesc(Y)
% colorbar
% colormap('jet')
% return

% adjust for cutoff in data
row2=500;
col2=row2;

x=zeros(row2,col2); % x coordinate of data map in pixels
y=zeros(row2,col2); % y coordinate of data map in pixels

yo=251;
xo=yo;

for i=1:col2
    x(:,i)=(-xo+i).*pix_size;  % counted from the left
end
for i=1:row2
    y(i,:)=(row2-yo+1-i)./cosd(ia).*pix_size; % counted from the bottom
end

mask=vy./vy;
mask(isnan(mask)==1)=0;

X_list=reshape(X,row*col,1);
Y_list=reshape(Y,row*col,1);
I_list=reshape(I,row*col,1);
vy_list=reshape(vy,row*col,1);
mask_list=reshape(mask,row*col,1);

X_mask_list=X_list;
Y_mask_list=Y_list;

X_list(isnan(vy_list)==1,:)=[];
Y_list(isnan(vy_list)==1,:)=[];
I_list(isnan(vy_list)==1,:)=[];
vy_list(isnan(vy_list)==1,:)=[];

% "interpolation"
F_I=TriScatteredInterp(X_list,Y_list,I_list,'natural');
F_vy=TriScatteredInterp(X_list,Y_list,vy_list,'natural');
F_mask=TriScatteredInterp(X_mask_list,Y_mask_list,mask_list,'natural');

mask=F_mask(x,y);
mask(mask<0.5)=nan;
mask(mask>=0.5)=1;

I_rot=F_I(x,y).*mask;
vy_rot=F_vy(x,y).*mask;

% save step 2 data
% on home pc
save('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\matlab_data\rotate_data.mat',....
    'I_rot','vy_rot','x','y');
% on mlc mac w/ usb
% save('/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/matlab_data/rotate_data.mat',....
%     'xo','yo','ia','pa','vsys','r25');


% write fits file for rotated data
% onto home pc
fitswrite(flipud(I_rot),'C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\fits_data\I_rot.fits');
fitswrite(flipud(vy_rot),'C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\fits_data\vy_rot.fits');
% on mlc mac w/ usb
% fitswrite(flipud(vy_rot),'/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/fits_data/vy_rot.fits');
% fitswrite(flipud(I_rot),'/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/fits_data/I_rot.fits');

% completed
goodbye=1

