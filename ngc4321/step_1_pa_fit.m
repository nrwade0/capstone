%/ 
% 12.12.2017
% Step 1: pa_fit NGC 4321
% fits for a more precise position angle based on velocity data
% outputs: figure 1... wide range (20) pa fit
% figure 2... small range (1) pa fit
% pa = 147.7
% vsys = 1559.440 
%/

clc
clear
format long g

% NGC 4321 initial values
xo=328.5;
yo=206.5;
ia=31.7;
pa=153;
r25=222.4;

% pix size found in DS9
pix_size=0.000223611*3600;

% velocity h_alpha data
% home pc path
data_path=fitsread('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\fits_data\n4321_v.fits'); 
% mlc lab mac on usb path
% data_path=fitsread('/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/fits_data/n4321_v.fits'); 

% correct data by flipping up/down
data=flipud(data_path);

% check for data
% figure(1)
% clf
% imagesc(data_path)
% colorbar
% colormap('jet')
% return

[row, col]=size(data);

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

x_list(data_list==0)=[];
y_list(data_list==0)=[];
data_list(data_list==0)=[];

x_list(isnan(data_list))=[];
y_list(isnan(data_list))=[];
data_list(isnan(data_list))=[];

% ring size
dr=5;

% ring dispersion from 80-100
r_ring=[80:dr:195]';
r_ring(:,2)=[85:dr:200]';

% number of rings found by number of rows in r_ring
n_ring=length(r_ring);

% increment pa's +- 20 and increment by 1
pa=[pa-20:1:pa+20]';

n=length(pa);
ssr=zeros(n,1);

for i=1:n
    i
    for j=1:n_ring
        
    X=(-(x_list-xo).*sind(pa(i,1))+(y_list-yo).*cosd(pa(i,1))).*pix_size; 
    Y=(((x_list-xo).*cosd(pa(i,1))+(y_list-yo).*sind(pa(i,1))))./cosd(ia).*pix_size; 
    R=sqrt(X.^2+Y.^2);

    d=data_list;
    
    d(R<=r_ring(j,1))=[];
    X(R<=r_ring(j,1))=[];
    Y(R<=r_ring(j,1))=[];
    R(R<=r_ring(j,1))=[];
    
        d(R>r_ring(j,2))=[];
    X(R>r_ring(j,2))=[];
    Y(R>r_ring(j,2))=[];
    R(R>r_ring(j,2))=[];
    
    W=abs(X./R);
    
    D=d.*W;
    
    g=X./R;
    
    G=ones(length(d),1).*W;
    
    G(:,2)=g.*W.*sind(ia);
    
    m=G\D;
    
    ssr(i,1)=ssr(i,1)+sum((D-G*m).^2);
    
    end
end

% plot of pa vs ssr
% figure(1)
% clf
% plot(pa,ssr,'b-')

pa_results=pa;
ssr_results=ssr;

% error: index exceeds matrix dimensions...
% all ssr's are = 0, so first run-thru (when pa=133)
% sets pa = 133 but it will run through the if condition since the first
% and second iteration are the same and try to set pa = 134
for i=1:n
    if ssr(i,1)==min(ssr)
        pa=pa(i,1)
    end
end

% go thru pa's +-1 by 0.1 increments
pa=[pa-1:0.1:pa+1]';

n=length(pa);
ssr=zeros(n,1);

for i=1:n
    i
    for j=1:n_ring
        
    X=(-(x_list-xo).*sind(pa(i,1))+(y_list-yo).*cosd(pa(i,1))).*pix_size; 
    Y=(((x_list-xo).*cosd(pa(i,1))+(y_list-yo).*sind(pa(i,1))))./cosd(ia).*pix_size; 
    R=sqrt(X.^2+Y.^2);

    d=data_list;
    
    d(R<=r_ring(j,1))=[];
    X(R<=r_ring(j,1))=[];
    Y(R<=r_ring(j,1))=[];
    R(R<=r_ring(j,1))=[];
    
    d(R>r_ring(j,2))=[];
    X(R>r_ring(j,2))=[];
    Y(R>r_ring(j,2))=[];
    R(R>r_ring(j,2))=[];
    
    W=abs(X./R);
    
    D=d.*W;
    
    g=X./R;
    
    G=ones(length(d),1).*W.*sind(ia);
    
    G(:,2)=g.*W;
    
    m=G\D;
    
    ssr(i,1)=ssr(i,1)+sum((D-G*m).^2);
    
    end
    
end

% second plot of pa vs ssr with smaller pa range
figure(2)
clf
plot(pa,ssr,'b-')
xlabel('Position Angle (\circ)');
ylabel('SSR (x10^{6})');

save('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\figs\pa_fit_fig.mat');

for i=1:n
    if ssr(i,1)==min(ssr);
        pa=pa(i,1)
    end
end

% calculation of vsys
X=(-(x_list-xo).*sind(pa)+(y_list-yo).*cosd(pa)); 
Y=(((x_list-xo).*cosd(pa)+(y_list-yo).*sind(pa)))./cosd(ia); 
R=sqrt(X.^2+Y.^2);

d=data_list;
    
d(R==0)=[];
X(R==0)=[];
Y(R==0)=[];
R(R==0)=[];
    
W=abs(X./R);
    
D=d.*W;
    
g=X./R;
    
G=ones(length(d),1).*W;
    
G(:,2)=g.*W.*sind(ia);
    
m=G\D;

% output vsys calc
vsys=m(1,1)

% save step 1 data onto usb
% on home pc
 save('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\matlab_data\pa_fit.mat',....
    'xo','yo','ia','pa','vsys','r25','pa_results','ssr_results')
% on mlc mac w/ usb
% save('/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/matlab_data/pa_fit.mat',....
%    'xo','yo','ia','pa','vsys','r25')


    
    



