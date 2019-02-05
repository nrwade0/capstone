clear
clc
format long g

% initial values
% xo, yo will change depending on data
vsys=1560.778933;  % from previous JK methods
xo=328.5;  
yo=206.5;  
ia=31.7;  
% pa=150;
pa=153;

pix_size=0.000223611*3600;
% beam_max=10.22;
% beam_min=9.07;
nc=5;

% implement velocity and intensity h_alpha data and correct vsys and ia
data_path=flipud(fitsread('C:\Users\Nick\Documents\MATLAB\h_alpha\n4321\fits_data\n4321_v.fits'));
vy=(data_path-vsys)./sind(ia);
% data_path='C:\Users\Nick\Documents\MATLAB\h_alpha\n4321\data\n4321_v.fits'; 
% vy=(flipud(1e-3.*fitsread(sprintf(data_path)))-vsys)./sind(ia);

data_path=fitsread('C:\Users\Nick\Documents\MATLAB\h_alpha\n4321\fits_data\n4321_i.fits');  
I=flipud(data_path);
% I=flipud(fitsread(sprintf(data_path)));

[row,col]=size(I);

% x, y coordinate of data map in pixels
x=zeros(row,col);
y=zeros(row,col);

% fill in coordinates for data maps
for i=1:col
    x(:,i)=i;  % counted from the left
end
for i=1:row
    y(i,:)=row-i+1; % counted from the bottom
end

% list of x, y coordinates in galaxy plane
X=(-(x-xo).*sind(pa)+(y-yo).*cosd(pa))*pix_size;
Y=((-(x-xo).*cosd(pa)-(y-yo).*sind(pa)))./cosd(ia)*pix_size;
R=sqrt(X.^2+Y.^2);

% redreate new data map
% x,y coordinate of other data map in pixels
x=zeros(row,col); 
y=zeros(row,col); 

for i=1:col
    x(:,i)=(-xo+i).*pix_size;  % counted from the left
end
for i=1:row
    y(i,:)=(row-yo+1-i)./cosd(ia).*pix_size; % counted from the bottom
end

% what does mask mean/do?
mask=vy./vy;
mask(isnan(mask)==1)=0;

% error in vy_list... cannot reshape because num of elements are changing.
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

F_I=TriScatteredInterp(X_list,Y_list,I_list,'natural');
F_vy=TriScatteredInterp(X_list,Y_list,vy_list,'natural');
F_mask=TriScatteredInterp(X_mask_list,Y_mask_list,mask_list,'natural');

mask=F_mask(x,y);

mask(mask<=0.5)=nan;
mask(mask>0.5)=1;

I_rot=F_I(x,y).*mask;
vy_rot=F_vy(x,y).*mask;

save('C:\Users\Nick\Documents\MATLAB\h_alpha\n4321\matlab_data\rotate_fit.mat',....
    'I_rot','vy_rot','x','y');
  
fitswrite(flipud(I_rot),'C:\Users\Nick\Documents\MATLAB\h_alpha\n4321\fits_data\I_rot.fits');
fitswrite(flipud(vy_rot),'C:\Users\Nick\Documents\MATLAB\h_alpha\n4321\fits_data\vy_rot.fits');

goodbye=1



