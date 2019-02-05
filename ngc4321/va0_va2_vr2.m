clc
clear
format long g

% initial values for the h_alpha data
vsys=1560.778933; % caught from JK earlier
xo=328.5;  
yo=206.5;  
ia=31.7; 
pa=153+1;
% pa=150;
r25=222.4;

% pa of the bar - found in profile_s2_3_6.m => phi_mean
phi=57.4743; % with pa=153
% phi=54.6451; % with pa=150
% phi=128.1884; % with pa=158
% phi=127.2398; % with pa=148
% phi=117.2296; % with pa=159
% phi=121.6092; % with pa=154

% input h_alpha velocity data map and correct for vsys and ia
data_path=fitsread('C:\Users\Nick\Documents\MATLAB\h_alpha\n4321\fits_data\n4321_v.fits');
data=((flipud(data_path))-vsys)./sind(ia);

% output figure to confirm similarity with DS9 image
% figure(1)
% clf
% imagesc(data)
% colorbar
% return

[row, col]=size(data);

mask=nan.*ones(row,col);

% pix_size found in display header of n4321_v.fits
% multiply by 3600 to convert into arcseconds
pix_size=0.000223611*3600;
nc=25; % 5?

x=zeros(row,col); 
y=zeros(row,col); 

for i=1:col
    x(:,i)=i;  
end
for i=1:row
    y(i,:)=row-i+1; 
end

x_list=reshape(x,row*col,1);
y_list=reshape(y,row*col,1);
data_list=reshape(data,row*col,1);

x_list(isnan(data_list))=[];
y_list(isnan(data_list))=[];
data_list(isnan(data_list))=[];

X=(-(x_list-xo).*sind(pa)+(y_list-yo).*cosd(pa)).*pix_size; 
Y=-(((x_list-xo).*cosd(pa)+(y_list-yo).*sind(pa)))./cosd(ia).*pix_size; 
R=sqrt(X.^2+Y.^2);

% size of data rings
dr=5;
% dr=10;

n=floor(r25./dr)-1;

r_ring=[dr:dr:n*dr]';
r_ring(:,2)=[2*dr:dr:(n+1)*dr]';

mean_r_ring=mean(r_ring')';

[va0, va2, vr2, k, va0_se, va2_se, vr2_se]=deal(zeros(n,1));

n % output length of loop
for i=1:n
    
    i % iterator
    
    d=data_list;
    x=X;
    y=Y;
    r=R;
    
    d(r<=r_ring(i,1))=[];
    x(r<=r_ring(i,1))=[];
    y(r<=r_ring(i,1))=[];
    r(r<=r_ring(i,1))=[];
    
    d(r>r_ring(i,2))=[];
    x(r>r_ring(i,2))=[];
    y(r>r_ring(i,2))=[];
    r(r>r_ring(i,2))=[];
    theta=atan2(y,x).*180./pi;

    G=x./r;
    G(:,2)=-cosd(2.*(theta-phi)).*x./r;
    G(:,3)=-sind(2.*(theta-phi)).*y./r;
    
    m=G\d;
    
    va0(i,1)=m(1,1);
    va2(i,1)=m(2,1);
    vr2(i,1)=m(3,1);
    
    k(i,1)=length(d);
    
    [va0_jk, va2_jk, vr2_jk]=deal(zeros(k(i,1),1));
    
    for j=1:k(i,1)
        g_jk=G;
        data_jk=d;
        
        data_jk(j,:)=[];
        g_jk(j,:)=[];
        
        m_jk=g_jk\data_jk;
        
        va0_jk(j,1)=m_jk(1,1);
        va2_jk(j,1)=m_jk(2,1);
        vr2_jk(j,1)=m_jk(3,1);
    end
    
    va0_se(i,1)=sqrt((k(i,1)/nc-1)/(k(i,1)/nc)*....
        sum((va0_jk-mean(va0_jk)).^2));
    va2_se(i,1)=sqrt((k(i,1)/nc-1)/(k(i,1)/nc)*....
        sum((va2_jk-mean(va2_jk)).^2));
    vr2_se(i,1)=sqrt((k(i,1)/nc-1)/(k(i,1)/nc)*.....
        sum((vr2_jk-mean(vr2_jk)).^2));
end

va0_ci=tinv(0.975,k./nc-1).*va0_se;
va2_ci=tinv(0.975,k./nc-1).*va2_se;
vr2_ci=tinv(0.975,k./nc-1).*vr2_se;

r=mean(r_ring')';

wp=(va0-va2)./r;

% plot of va0 and their respective confidence intervals
figure(1)
clf
errorbar(r,va0, va0_ci,'b*')

% should cross 0 line where bar is located
figure(2)
clf
errorbar(r,va2, va2_ci,'bo')
hold on
plot([0 r25],[0 0],'k-')
% plot([220 220],[-40 10],'k-')
hold off

% should cross 0 line where bar is located
figure(3)
clf
errorbar(r,vr2,vr2_ci,'ro')
hold on
plot([0 r25],[0 0],'k-')
% plot([220 220],[-40 10],'k-')
hold off

% model of radial flow and pattern speed
% should be 'bump' in flow around the location of the bar
figure(4)
clf
plot(r,va0./r,'k-')
hold on
plot(r,wp,'r-')
hold off
axis([0 250 1 5])





