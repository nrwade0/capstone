clear
clc

% open on mlc mac w/ usb
% load('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\figs\reg_fit_all_o.mat');
% open on home pc
load('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\figs\step4_figs\step4_fig_n11.mat');

% axis ticks
xtick_label1=[0 20 40 60 80 100 120 140 160 180];
xtick1=[0 20 40 60 80 100 120 140 160 180];
ytick_label1=[0 1 2 3 4 5 6 7];
ytick1=[0 1 2 3 4 5 6 7];

fontsize1=10;
fontsize2=9;
fontsize3=9;

L1=0.5;

figure(1)
clf

% x label
subplot('position',[0.199 0.075 0.01 0.01])
set(gca,'FontName','Arial','FontSize',fontsize1,...
    'TickLength', [0.02 0.02],'XMinorTick','off',...
    'YMinorTick','off','xtick',[],'xticklabel',[],...
    'ytick',[],'yticklabel',[],'XColor','white','YColor','white') 
xlabel('$r$ \sf{(arcsec)}','Fontsize',fontsize1,'Interpreter','Latex','Color','black')

% y label units only
% subplot('position',[0.076 0.174 0.001 0.001])
% set(gca,'FontName','Arial','FontSize',fontsize1,...
%     'TickLength', [0.02 0.02],'XMinorTick','off',...
%     'YMinorTick','off','xtick',[],'xticklabel',[],...
%     'ytick',[],'yticklabel',[],'XColor','white','YColor','white') 
% ylabel('\sf{(km s$^{-1}$ arcsec$^{-1}$)}','Fontsize',fontsize1,'Interpreter','Latex','Color','black')

% y label words
subplot('position',[0.05 0.17 0.001 0.001])
set(gca,'FontName','Arial','FontSize',fontsize1,...
    'TickLength', [0.02 0.02],'XMinorTick','off',...
    'YMinorTick','off','xtick',[],'xticklabel',[],...
    'ytick',[],'yticklabel',[],'XColor','white','YColor','white') 
ylabel('\sf{Angular Frequency} \sf{(km s$^{-1}$ arcsec$^{-1}$)}','Fontsize',fontsize1,'Interpreter','Latex','Color','black')

[r_terpw, r_terp2m, r_terp2p, r_terp4p, r_terp4m]=deal(r_terp);

% resonances
l2m=w-K/2;
l2p=w+K/2;
l4m=w-K/4;
l4p=w+K/4;

% cutoff for resonances
cut = 7.5;

% set cutoff for resonances, etc
r_terpw(w>cut)=[];
w(w>cut)=[];

r_terp2m(l2m>cut)=[];
l2m(l2m>cut)=[];

r_terp2p(l2p>cut)=[];
l2p(l2p>cut)=[];

r_terp4m(l4m>cut)=[];
l4m(l4m>cut)=[];

r_terp4p(l4p>cut)=[];
l4p(l4p>cut)=[];

% wtf are these for
ye=1.8;
yb=0;

% subplot('position',[0.1 0.1 0.22 0.355*(ye-yb)/3.7]);
subplot('position',[0.1 0.1 0.85 0.85]);

% add confidence intervals 
for i=1:n_ring-3
    top=wp(i,1)+wp_ci(i,1);
    bott=wp(i,1)-wp_ci(i,1);
    fill([r_ring(i,1) r_ring(i,2) r_ring(i,2) r_ring(i,1)],....
        [top top bott bott],[1,0.80,0.80],....
        'edgecolor',[1,0.8,0.8])
    if i==1
        hold on
    end
end

% plot resonances, etc
plot(r_terpw,w,':k','LineWidth',1)
plot(r_terp2m,l2m,'b--','LineWidth',1)
plot(r_terp2p,l2p,'b--','LineWidth',1)
plot(r_terp4p,l4p,'g-.','LineWidth',1)
plot(r_terp4m,l4m,'g-.','LineWidth',1)

% plot omega p
for i=1:n_ring-3
    % plot omega p in the middle of ci
%     plot([r_ring(i,1) r_ring(i,2)],wp(i,1).*ones(2,1),'-','Color',[0.7 0 0],'LineWidth',1.5)
    
    % plot ci bands w/ omega p on top and bottom of ci
    plot([r_ring(i,1) r_ring(i,2)],wp(i,1).*ones(2,1)+wp_ci(i,1),'-','Color',[0.7 0 0],'LineWidth',1.5)
    plot([r_ring(i,1) r_ring(i,2)],wp(i,1).*ones(2,1)-wp_ci(i,1),'-','Color',[0.7 0 0],'LineWidth',1.5)
end
hold off

axis([10 170 yb ye])
% first 2: x-axis; last 2: y-axis
axis([10 170 0 7])

set(gca,'FontName','Arial','FontSize',fontsize2,...
    'TickLength', [0.02 0.02],'XMinorTick','on',...
    'YMinorTick','on','xtick',xtick1,'xticklabel',xtick_label1,...
    'ytick',ytick1,'yticklabel',ytick_label1)

% adds small figure "(b)" in bottom left of figure
% text(7, 0.15,'($b$)','Fontsize',fontsize3,'Interpreter', 'Latex')

a='C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\figs\finished\wp_all_n11.eps';
set(figure(1), 'color', 'white')
print(figure(1),'-r600', '-depsc', a)
fix_lines(a)

goodbye=1
