%/
% Step 4: reg_fit
% input: step 1-3 data, pix size
% output: many tests
% final graph - muy importante
%/

clear
clc
format

% pix size from DS9
pix_size=0.000223611*3600;

% load data from steps 1 and 3
% on home pc
load('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\matlab_data\pa_fit.mat')
load('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\matlab_data\vel_fit.mat')
% on mlc mac w/ usb
% load('/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/matlab_data/pa_fit.mat')
% load('/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/matlab_data/vel_fit.mat')

% reset yo to adjust for rotation
yo=251;
xo=251;

r_terp=[r_ring(1,1):1:r_ring(length(r_ring),2)]';
v_terp=interp1(mean(r_ring')',v,r_terp,'pchip');

w=v_terp./r_terp;

diff_wsqr=-v_terp.^2./r_terp.^3;

K=sqrt(r_terp.*diff_wsqr+4.*w.^2);

% load step 2 data
% on home pc
I=flipud(fitsread('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\fits_data\I_rot.fits'));
vy=flipud(fitsread('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\fits_data\vy_rot.fits'));
% on mlc mac w/ usb
% I=flipud(fitsread('/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/fits_data/I_rot.fits'));
% vy=flipud(fitsread('/Volumes/Untitled/MATLAB_prefail/h_alpha/n4321/fits_data/vy_rot.fits'));

% % make sure step 2 data looks ok
% figure(1)
% clf
% imagesc(vy)
% colormap('jet')
% colorbar
% return
% 
% figure(2)
% clf
% imagesc(I)
% colormap('jet')
% colorbar
% return

[row,col]=size(I);

[x, y]=deal(zeros(row,col));

for i=1:col
    x(:,i)=(i-xo).*pix_size;
end

for i=1:row
    y(i,:)=(row-yo-i+1)/cosd(ia).*pix_size;
end

Ix=I.*x;
Ivy=I.*vy;

Ix(isnan(Ix))=0;
Ivy(isnan(Ivy))=0;

r=sqrt(x.^2+y.^2);

% edit
kill_ring=10;

r_ring(length(r_ring)-kill_ring:length(r_ring),:)=[];
r_ring=[0 r_ring(1,1); r_ring; r_ring(length(r_ring),2) 1e3];
n_ring=length(r_ring);
r_ring_mean=mean(r_ring')';

% edit edit edit edit edit edit edit edit edit edit edit edit edit edit 
n=0;
% n=5;
% n=11;

if n~=0
    
lim1=r_ring(n,2);

for i=1:n
    r_ring(1,:)=[];
end

n_ring=length(r_ring);

else
    lim1=0;
end

% edit?
lim2=2000;

for i=1:n_ring
    % i
    
    Ix_ring=zeros(row,col);
    Ivy_ring=zeros(row,col);
        
    for j=1:row
       for k=1:col
           if abs(y(j,1))>=lim1 
                if abs(y(j,1))<lim2;
                    if r(j,k)>r_ring(i,1) && r(j,k)<=r_ring(i,2)
                    Ix_ring(j,k)=Ix(j,k);
                    Ivy_ring(j,k)=Ivy(j,k);
                    end
               end
           end
        end
    end
        
    % watch rings form 
    % shows what?
%     figure(1)
%     clf
%     imagesc(Ix_ring)
%     pause(1)
    
    G(:,i)=sum(Ix_ring,2);
    d_map(:,i)=sum(Ivy_ring,2);

    % watch values of G progress thru rings
    % shows what?
%     figure(2)
%     clf
%     imagesc(G)
%     pause(1)
end

d=sum(d_map,2);

d_sum=d;
G_test=sum(G,2);

% % test for what
% figure(1)
% clf
% imagesc(log(abs(G)))
% colormap('jet')
% colorbar

% edit
n=249;
g=zeros(n,length(r_ring));
g_test=zeros(n,length(r_ring));
d=zeros(length(r_ring),1);
    
for i=1:n
    g(i,:)=G(i,:)+G(row-i,:);
    d(i,:)=d_sum(i,:)+d_sum(row-i,:);
    g_test(i,:)=G_test(i,:)+G_test(row-i,:);
end

G=g;
G_test=g_test;

% % test for what
% figure(2)
% clf
% imagesc(log(abs(G)))
% colormap('jet')
% colorbar

G_test=sum(G,2);
k=0;

for j=1:length(G_test)
    o=j-k;
    if G_test(o,1)==0
        k=k+1;
        G_test(o,:)=[];
        G(o,:)=[];
        d(o,:)=[];
    end
end

L=get_l(n_ring,1);

[UU,sm,XX]=cgsvd(G,L);

figure(1)
clf
[reg_corner,rho,eta,reg_param] = l_curve(UU,sm,d,'Tikh');

% tikhonov regularization
% reg_corner is suggested based on model but
% 3e5 is our set value for a
a=reg_corner
a=3e5;

[wp,rho,eta]=tikhonov(UU,sm,XX,d,a);

n=length(G(:,1));

wp_jk=zeros(n_ring,n);

for i=1:n
    i
    G_jk=G;
    d_jk=d;

    G_jk(i,:)=[];
    d_jk(i,:)=[];
    
    [UU,sm,XX]=cgsvd(G_jk,L);
    
    [reg_corner_jk,rho_jk,eta_jk,reg_param] = l_curve(UU,sm,d_jk,'Tikh');

%     a=reg_corner;
    
    [wp_jk(:,i),rho,eta]=tikhonov(UU,sm,XX,d_jk,a) ;

end

[wp_ci, wp_se]=deal(zeros(n_ring,1));

% calculate confidence intervals for omega p
for i=1:n_ring
    i
    wp_se(i,1)=sqrt(sum((mean(wp_jk(i,:))-wp_jk(i,:)).^2)*(n-1)/(n));
    wp_ci(i,1)=tinv(0.975,n-1).*wp_se(i,1);

end

% the most important figure in the project
figure(2)
clf
% plot omega
plot(r_terp,w,':k')
hold on
% plot inner lindblad resonance
plot(r_terp,w+K/4,'g-.')
plot(r_terp,w-K/4,'g-.')
% plot outer lindblad resonance
plot(r_terp,w+K/2,'b:')
plot(r_terp,w-K/2,'b:')

for i=1:n_ring-1
%   plot([r_ring(i,1) r_ring(i,2)],(wp(i,1)+wp_ci(i,1)).*ones(2,1),'r-')
%   plot([r_ring(i,1) r_ring(i,2)],(wp(i,1)-wp_ci(i,1)).*ones(2,1),'r-')
    plot([r_ring(i,1) r_ring(i,2)],(wp(i,1)).*ones(2,1),'r-')
end
hold off

% add labels and axis dimensions
xlabel('r (arcsec)');
ylabel('Angular Frequency (km s^{-1} arsec^{-1})');
axis([0 190 0 7])

% save step 4 fig
save('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\figs\step4_figs\step4_fig_n11.mat',...
    'wp','wp_ci','r_ring','r_terp','w','K','n_ring')

% add legend in top right
% legend([\omega + \kappa / 4, Fit, Data, Model, ], 'Data (\mu \pm \sigma)', ...
%    'Fit (\it{C x^3})', 'Validation Data', 'Model (\it{C x^3})', '95% CI', ...
%    'location', 'NorthEast');

% calculate mean of bar speedwhen
% n=3
% b=mean(wp(1:9,1))
% n=5
% b=mean(wp(1:7,1))

a=mean(wp(1:11));
se=std(wp(1:11)/sqrt(11));
ci=tinv(0.975,10)*se;



