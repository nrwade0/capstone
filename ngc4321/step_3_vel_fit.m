%/
% Step 3: vel_fit
% inputs: pa_fit data, velocity map (not rotated)
% outputs: velocity vs v confidence interval
% angular frequency graphs
% fig 1 = ln curve, fig 2 = 1/x curve
%/

clear
clc
format long g

% load step 1 data
% on home pc
load('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\matlab_data\pa_fit.mat');
% onto mlc mac w/ usb
% load('/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/matlab_data/pa_fit.mat');

% input velocity map
% on home pc
data_path=fitsread('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\fits_data\n4321_v.fits'); 
% on mlc mac w/ usb
% data_path=fitsread('/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/fits_data/n4321_v.fits');

% data corrections: flips and converts to km/s
data=flipud(data_path);
% set data < 0 to nan
data(data==0)=nan;
% correct for inclination angle and systemic velocity
data=(data-vsys)./sind(ia);

% creates data map x,y via row,col
[row, col]=size(data);
[x, y]=deal(zeros(row,col));

for i=1:col 
    x(:,i)=i;
end

for i=1:row
    y(i,:)=(row-i+1);
end

data=reshape(data,row*col,1);
x=reshape(x,row*col,1);
y=reshape(y,row*col,1);

x(isnan(data))=[];
y(isnan(data))=[];
data(isnan(data))=[];

% data ring size (used 10 for n4736)
dr=5;
% dr=10;

n_ring=floor(r25/dr)-1;

r_ring=[dr:dr:n_ring*dr]';
r_ring(:,2)=[2*dr:dr:(n_ring+1)*dr]';

% number of rings is the amount of rings in r_ring matrix
n_ring=length(r_ring);

% make galaxy coordinates
X=-(x-xo).*sind(pa)+(y-yo).*cosd(pa);
Y=-((x-xo).*cosd(pa)+(y-yo).*sind(pa))./cosd(ia);
R=sqrt(X.^2+Y.^2);

[v, v_ci]=deal(zeros(n_ring,1));

for i=1:n_ring
    
    i
    W=abs(X./R);
    
    d=data.*W;
    cos_theta=X./R.*W;
    r=R;
    
    d(r<=r_ring(i,1),:)=[];
    cos_theta(r<=r_ring(i,1),:)=[];
    r(r<=r_ring(i,1),:)=[];
    
    d(r>r_ring(i,2),:)=[];
    cos_theta(r>r_ring(i,2),:)=[];
    r(r>r_ring(i,2),:)=[];
    
    v(i,1)=cos_theta\d;
    
    v_jk=zeros(length(d),1);
    
    for j=1:length(d)
        
        d_jk=d;
        cos_theta_jk=cos_theta;
        
        d_jk(j,:)=[];
        cos_theta_jk(j,:)=[];
        
        v_jk(j,1)=cos_theta_jk\d_jk;
    end
        
    v_ci(i,1)=tinv(0.975,length(d)-1)*sqrt(sum((mean(v_jk)-v_jk).^2)*(length(d)-1)/length(d));
    
end

% velocity makes ln curve
% what do they show?
figure(1)
clf
errorbar(mean(r_ring')',v,v_ci,'ro')

% angular frequency makes 1/x shape
% what do they show?
% figure(2)
% clf
% plot(mean(r_ring')',v./mean(r_ring')','b*')

% save step 3 data
% onto home pc
save('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\matlab_data\vel_fit.mat',...
    'r_ring','v','v_ci')
% on mlc mac w/ usb
% save('/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/matlab_data/vel_fit.mat',...
%     'r_ring','v','v_ci')


