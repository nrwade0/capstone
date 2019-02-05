clear
clc
format

% input rotated 3.6 micrometer fits data of n4321
data=flipud(fitsread('C:\Users\Nick\Documents\MATLAB\h_alpha\n4321\fits_data\n4321_rotated_3_6.fits'));

[row, col]=size(data);

% intital values given in parameters.txt
xo=721.851;
yo=1205.553;  
ia=31.7; 
pa=150;
% pa=153;
r25=222.4;

% for n4321_rotated_3_6.fits pix_size=0.75
pix_size=0.75;
nc=25;

% size of rings
dr=5;
%dr=10;

n_ring=floor(r25/dr)-1;
% n_ring=33;

r_ring=[dr:dr:(n_ring)*dr]';
r_ring(:,2)=[2*dr:dr:(n_ring+1)*dr]';
r_ring_mean=mean(r_ring')';

x=zeros(row,col); 
y=zeros(row,col);

for i=1:col
    x(:,i)=i;  
end
for i=1:row
    y(i,:)=row-i+1; 
end

X=(-(x-xo).*sind(pa)+(y-yo).*cosd(pa)).*pix_size; 
Y=-(((x-xo).*cosd(pa)+(y-yo).*sind(pa)))./cosd(ia).*pix_size; 
R=sqrt(X.^2+Y.^2);

Theta=atan2(Y,X).*180./pi;

X_list=reshape(X,row*col,1);
Y_list=reshape(Y,row*col,1);
Theta_list=reshape(Theta,row*col,1);
data_list=reshape(data,row*col,1);

R_list=sqrt(X_list.^2+Y_list.^2); 

Theta_list(isnan(data_list)==1)=[];
R_list(isnan(data_list)==1)=[];
data_list(isnan(data_list)==1)=[];

Theta_list(data_list==0)=[];
R_list(data_list==0)=[];
data_list(data_list==0)=[];

R_list(Theta_list<0)=[];
data_list(Theta_list<0)=[];
Theta_list(Theta_list<0)=[];

[I0, I2c, I2s, I0_se, I2c_se, I2s_se, k]=deal(zeros(n_ring,1));

[D, T, M]=deal(cell(n_ring,1));

for i=1:n_ring
    i  % iterator
    
    Theta=Theta_list;
    data=data_list;
    R=R_list;
    
    Theta(R<=r_ring(i,1))=[];
    data(R<=r_ring(i,1))=[];
    R(R<=r_ring(i,1))=[];
    
    Theta(R>r_ring(i,2))=[];
    data(R>r_ring(i,2))=[];
    R(R>r_ring(i,2))=[];
    
    g=ones(length(data),1);
    g(:,2)=cosd(2*Theta);
    g(:,3)=sind(2*Theta);

    m=g\data;
    
    I0(i,1)=m(1,1);
    I2c(i,1)=m(2,1);
    I2s(i,1)=m(3,1);
    
    k(i,1)=length(data);
    
    [I0_jk, I2c_jk, I2s_jk]=deal(zeros(k(i,1),1));
    
    for j=1:k(i,1)
        g_jk=g;
        data_jk=data;
        
        data_jk(j,:)=[];
        g_jk(j,:)=[];
        
        m_jk=g_jk\data_jk;
        
        I0_jk(j,1)=m_jk(1,1);
        I2c_jk(j,1)=m_jk(2,1);
        I2s_jk(j,1)=m_jk(3,1);
    end

    I0_se(i,1)=sqrt((k(i,1)/nc-1)/(k(i,1)/nc)*....
        sum((I0_jk-mean(I0_jk)).^2));
    I2c_se(i,1)=sqrt((k(i,1)/nc-1)/(k(i,1)/nc)*....
        sum((I2c_jk-mean(I2c_jk)).^2));
    I2s_se(i,1)=sqrt((k(i,1)/nc-1)/(k(i,1)/nc)*.....
        sum((I2s_jk-mean(I2s_jk)).^2));
    
    D{i,1}=data;
    T{i,1}=Theta;
    M{i,1}=m;    
end

% end results

phi_3_6=0.5.*atan2(I2s,I2c).*180/pi;

phi_3_6_se=0.5.*sqrt((I2s_se./I2c).^2+ (I2c_se.*I2s./I2c.^2).^2)./....
    (1+(I2s./I2c).^2).*180./pi;

phi_3_6_ci=tinv(0.975,k./nc-3).*phi_3_6_se;

% comment this out the first time and edit as you see fit to make the
% result look continuous
for i=1:n_ring
    if phi_3_6(i,1)<0;
        phi_3_6(i,1)=phi_3_6(i,1)+180;
    end
end

% calculate weighted mean of phi and output
g=1./phi_3_6_ci(1:13,1);
d=phi_3_6(1:13,1)./phi_3_6_ci(1:13,1);
phi_mean=g\d

% confidence interval of the mean of phi and output
phi_mean_ci=(std(phi_3_6(1:13,1)))./sqrt(13)

% large distortion in the center of the galaxy
% flat part is the bar
% 'drop off' is where the bar ends
figure(1)
clf
errorbar(mean(r_ring')',phi_3_6,phi_3_6_ci,'b*')
xlabel('r (arcsec)','Interpreter','LaTex')
ylabel('$\theta$ (degrees)','Interpreter','LaTex')
% axis([0 350 30 120]) % you'll want to comment this the first time and edit it

theta_plot=[0:1:179]';

% for i=1:n_ring
%     figure(i+1)
%     clf
%     plot(T{i,1},D{i,1},'b*')
%     hold on
%     m=M{i,1};
%     plot(theta_plot,m(1,1)+m(2,1).*cosd(2.*theta_plot)+m(3,1).*sind(2.*theta_plot),'r-')
%     hold off
% end
    

save('C:\Users\Nick\Documents\MATLAB\h_alpha\n4321\matlab_data\profile_s2_3_6_data.mat',....
    'phi_3_6','phi_3_6_ci','r_ring','phi_3_6_se','T','D','M','theta_plot')

return

load('C:\Users\Nick\Documents\MATLAB\h_alpha\n4321\matlab_data\profile_s2_3_6_data.mat');

figure(1)
clf
errorbar(mean(r_ring')',phi_3_6,phi_3_6_ci,'b*')
xlabel('r (arcsec)','Interpreter','LaTex')
ylabel('$\theta$ (degrees)','Interpreter','LaTex')
axis([0 350 30 110]) 



