clear
clc

% load given data
load('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\matlab_data\ellipse_fit.mat');

r_terp=[r_ring(1,1):1:r_ring(length(r_ring),2)]';
va0_terp=interp1(mean(r_ring')',va0,r_terp,'pchip');
w=va0_terp./r_terp;
diff_wsqr=-va0_terp.^2./r_terp.^3;
K=sqrt(r_terp.*diff_wsqr+4.*(w).^2);

[r_terpw, r_terp2m, r_terp2p, r_terp4p, r_terp4m]=deal(r_terp);

l2m=w-K/2;
l2p=w+K/2;
l4m=w-K/4;
l4p=w+K/4;

% cuts
cut1 = 25;
cut2 = 60;
cut3 = 6.4;

% apply cuts to figure lines
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

% x y labels and ticks
xtick_label1=[0 50 100 150];
xtick1=[0 50 100 150];
ytick_label1=[0 2 4 6];
ytick1=[0 2 4 6];

% font sizes (do not change)
f1=20;
f2=20;

figure(1)
clf

% add x label
subplot('position',[0.57 0.07 0.01 0.01])
set(gca,'FontSize',f1,...
    'TickLength', [0.02 0.02],'XMinorTick','off',...
    'YMinorTick','off','xtick',[],'xticklabel',[],...
    'ytick',[],'yticklabel',[],'XColor','white','YColor','white') 
xlabel('$r$ (arcsec)','Fontsize',f1,'Interpreter', 'Latex','Color','black')

% add y label (title)
subplot('position',[0.07 0.5 0.01 0.1])
set(gca,'FontSize',f1,...
    'TickLength', [0.02 0.02],'XMinorTick','off',...
    'YMinorTick','off','xtick',[],'xticklabel',[],...
    'ytick',[],'yticklabel',[],'XColor','white','YColor','white') 
ylabel('Angular Frequency','Fontsize',f1,'Interpreter', 'Latex','Color','black')

% add y label (units) underneath
subplot('position',[0.12 0.5 0.01 0.1])
set(gca,'FontSize',f1,...
    'TickLength', [0.02 0.02],'XMinorTick','off',...
    'YMinorTick','off','xtick',[],'xticklabel',[],...
    'ytick',[],'yticklabel',[],'XColor','white','YColor','white') 
ylabel('(km s$^{-1}$ arcsec$^{-1}$)','Fontsize',f1,'Interpreter', 'Latex','Color','black')

subplot('position',[0.16 0.13 0.83 0.86])

r_ring=r_ring(5:11,:);

for i=1:length(wp)
    
 top=wp(i,1)+wp_ci(i,1);
 bott=wp(i,1)-wp_ci(i,1);
 fill([r_ring(i,1) r_ring(i,2) r_ring(i,2) r_ring(i,1)],....
        [top top bott bott],[1,0.80,0.80],....
        'edgecolor',[1,0.8,0.8])
    
    if i==1
        hold on
    end
end

% plot figure lines
plot(r_terpw,w,':k','LineWidth',1)
plot(r_terp2m,l2m,'b--','LineWidth',1)
plot(r_terp2p,l2p,'b--','LineWidth',1)
plot(r_terp4p,l4p,'g-.','LineWidth',1)
plot(r_terp4m,l4m,'g-.','LineWidth',1)

top=3.17+0.04;
bott=3.17-0.04;
fill([r_ring(1,1) r_ring(length(wp),2) r_ring(length(wp),2) r_ring(1,1) ],....
    [top top bott bott],[0.8,0.8,0],....
    'edgecolor',[0.8,0.8,0])
    
for i=1:length(wp)
    plot([r_ring(i,1) r_ring(i,2)],(wp(i,1)+wp_ci(i,1)).*ones(2,1),'-','Color',[0.7 0 0],'LineWidth',1)
    plot([r_ring(i,1) r_ring(i,2)],(wp(i,1)-wp_ci(i,1)).*ones(2,1),'-','Color',[0.7 0 0],'LineWidth',1)
end

hold off

% match axes with fig_wp_vel
axis([-5 170 -0.2 6.9])

set(gca,'FontSize',f2,...
    'TickLength', [0.02 0.02],'XMinorTick','on',...
    'YMinorTick','on','xtick',xtick1,'xticklabel',xtick_label1,...
    'ytick',ytick1,'yticklabel',ytick_label1)

% save as eps file in fig finished folder
a='C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\figs\finished\wp_disc.eps';
set(figure(1), 'color', 'white')
print(figure(1),'-r600', '-depsc', a)


