clear
clc

load('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\matlab_data\pa_fit.mat');

% font sizes (do not change - 20)
f1=20;
f2=20;

% axis ticks
xtick_label1=[140 150 160 170];
xtick1=[140 150 160 170];
ytick_label1=[5 10 15 20 25];
ytick1=[5 10 15 20 25];

% slim down x10^6
ssr_results=ssr_results./10^6;

% x label (SSR units: x10^6 km^2 s^-2)
subplot('position',[0.4 0.7 0.2 0.01])
set(gca,'FontSize',f1,...
    'TickLength', [0.02 0.02],'XMinorTick','off',...
    'YMinorTick','off','xtick',[],'xticklabel',[],...
    'ytick',[],'yticklabel',[],'XColor','white','YColor','white') 
xlabel('Position Angle (degrees)','Fontsize',f1,'Interpreter', 'Latex','Color','black')

% y label
subplot('position',[0.2 0.5 0.01 0.1])
set(gca,'FontSize',f1,...
    'TickLength', [0.02 0.02],'XMinorTick','off',...
    'YMinorTick','on','xtick',[],'xticklabel',[],...
    'ytick',[],'yticklabel',[],'XColor','white','YColor','white') 
ylabel('SSR ($\times10^{6}$ km$^{-2}$ s$^{-2}$)','Fontsize',f1,'Interpreter', 'Latex','Color','black')

xlabel('Position Angle (degrees)','Fontsize',f1,'Interpreter', 'Latex','Color','black');
ylabel('SSR ($\times10^{6}$ km$^{-2}$ s$^{-2}$)','Fontsize',f1,'Interpreter', 'Latex','Color','black');

figure(1)
clf

% create drawing board for results plot
subplot('position',[0.15 0.15 0.81 0.81])

% plot pa_fit results
plot(pa_results,ssr_results,'r-');

% x-left x-right y-bottom y-top
axis([132 174 4 29])

set(gca,'FontSize',f2,...
    'TickLength', [0.02 0.02],'XMinorTick','on',...
    'YMinorTick','on','xtick',xtick1,'xticklabel',xtick_label1,...
    'ytick',ytick1,'yticklabel',ytick_label1)

% save data as eps file
a='C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\figs\finished\pa_fit.eps';
set(figure(1), 'color', 'white')
print(figure(1),'-r600', '-depsc', a)


