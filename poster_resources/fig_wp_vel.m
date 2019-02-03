clear
clc

% fix x location

load('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\matlab_data\vel_fit.mat');

% mismatched variable names
va0=v;

r_terp=[r_ring(1,1):1:r_ring(length(r_ring),2)]';
va0_terp=interp1(mean(r_ring')',va0,r_terp,'pchip');

w=va0_terp./r_terp;

diff_wsqr=-va0_terp.^2./r_terp.^3;

K=sqrt(r_terp.*diff_wsqr+4.*(w).^2);

[r_terpw, r_terp2m, r_terp2p, r_terp4p, r_terp4m]=deal(r_terp);

% resonances
l2m=w-K/2;
l2p=w+K/2;
l4m=w-K/4;
l4p=w+K/4;

% edit
% cut resonances to fit our dimensions of resonances
% cut1: left
% cut2: right
% cut3: up
cut1 = 0;
cut2 = 165;
cut3 = 6.75;

l2m(r_terpw>cut2)=[];
l2p(r_terpw>cut2)=[];
l4m(r_terpw>cut2)=[];
l4p(r_terpw>cut2)=[];
r_terp2m(r_terpw>cut2)=[];
r_terp2p(r_terpw>cut2)=[];
r_terp4m(r_terpw>cut2)=[];
r_terp4p(r_terpw>cut2)=[];
w(r_terpw>cut2)=[];
r_terpw(r_terpw>cut2)=[];

l2m(r_terpw<cut1)=[];
l2p(r_terpw<cut1)=[];
l4m(r_terpw<cut1)=[];
l4p(r_terpw<cut1)=[];
r_terp2m(r_terpw<cut1)=[];
r_terp2p(r_terpw<cut1)=[];
r_terp4m(r_terpw<cut1)=[];
r_terp4p(r_terpw<cut1)=[];
w(r_terpw<cut1)=[];
r_terpw(r_terpw<cut1)=[];

r_terp2m(l2m>cut3)=[];
l2m(l2m>cut3)=[];
r_terp2p(l2p>cut3)=[];
l2p(l2p>cut3)=[];
r_terp4m(l4m>cut3)=[];
l4m(l4m>cut3)=[];
r_terp4p(l4p>cut3)=[];
l4p(l4p>cut3)=[];
r_terpw(w>cut3)=[];
w(w>cut3)=[];

% mismatched variable names
W=w;

% change loaded data
load('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\figs\step4_figs\step4_fig_n11.mat');

% ticks
xtick_label1=[0 50 100 150];
xtick1=[0 50 100 150];
ytick_label1=[0 2 4 6];
ytick1=[0 2 4 6];

% font sizes (do not change)
f1=20;
f2=20;

figure(1)
clf

subplot('position',[0.57 0.07 0.01 0.01])
set(gca,'FontSize',f1,...
    'TickLength', [0.02 0.02],'XMinorTick','off',...
    'YMinorTick','off','xtick',[],'xticklabel',[],...
    'ytick',[],'yticklabel',[],'XColor','white','YColor','white') 
xlabel('$r$ (arcsec)','Fontsize',f1,'Interpreter', 'Latex','Color','black')

subplot('position',[0.07 0.5 0.01 0.1])
set(gca,'FontSize',f1,...
    'TickLength', [0.02 0.02],'XMinorTick','off',...
    'YMinorTick','off','xtick',[],'xticklabel',[],...
    'ytick',[],'yticklabel',[],'XColor','white','YColor','white') 
ylabel('Angular Frequency','Fontsize',f1,'Interpreter', 'Latex','Color','black')

subplot('position',[0.12 0.5 0.01 0.1])
set(gca,'FontSize',f1,...
    'TickLength', [0.02 0.02],'XMinorTick','off',...
    'YMinorTick','off','xtick',[],'xticklabel',[],...
    'ytick',[],'yticklabel',[],'XColor','white','YColor','white') 
ylabel('(km s$^{-1}$ arcsec$^{-1}$)','Fontsize',f1,'Interpreter', 'Latex','Color','black')

subplot('position',[0.16 0.13 0.83 0.86])

% edit
% r_ring=r_ring(1:34,:);

for i=1:length(wp)-1
     top=wp(i,1)+wp_ci(i,1);
     bott=wp(i,1)-wp_ci(i,1);
     fill([r_ring(i,1), r_ring(i,2), r_ring(i,2), r_ring(i,1)],....
        [top top bott bott],[1,0.80,0.80],....
        'edgecolor',[1,0.8,0.8])
    
    if i==1
        hold on
    end
end

% plot resonances and patten speed
plot(r_terpw,W,':k','LineWidth',1)
plot(r_terp2m,l2m,'b--','LineWidth',1)
plot(r_terp2p,l2p,'b--','LineWidth',1)
plot(r_terp4p,l4p,'g-.','LineWidth',1)
plot(r_terp4m,l4m,'g-.','LineWidth',1)
    
for i=1:length(wp)-1
    plot([r_ring(i,1) r_ring(i,2)],(wp(i,1)+wp_ci(i,1)).*ones(2,1),'-','Color',[0.7 0 0],'LineWidth',1)
    plot([r_ring(i,1) r_ring(i,2)],(wp(i,1)-wp_ci(i,1)).*ones(2,1),'-','Color',[0.7 0 0],'LineWidth',1)
end

hold off

% do not change axis
axis([-5 170 -0.2 6.9])

set(gca,'FontSize',f2,...
    'TickLength', [0.02 0.02],'XMinorTick','on',...
    'YMinorTick','on','xtick',xtick1,'xticklabel',xtick_label1,...
    'ytick',ytick1,'yticklabel',ytick_label1)

a='C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\figs\finished\n11_wp_all.eps';
set(figure(1), 'color', 'white')
print(figure(1),'-r600', '-depsc', a)



