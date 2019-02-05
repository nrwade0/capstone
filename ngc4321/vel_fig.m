clear
clc

load('C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\matlab_data\vel_fit.mat');

% font sizes (do not change)
f1=20;
f2=20;

% ticks
xtick_label1=[0 50 100 150];
xtick1=[0 50 100 150];
ytick_label1=[75 100 125 150 175 200 225 250];
ytick1=ytick_label1;

figure(1)
clf

subplot('position',[0.57 0.07 0.01 0.01])
set(gca,'FontSize',f1,...
    'TickLength', [0.02 0.02],'XMinorTick','off',...
    'YMinorTick','off','xtick',[],'xticklabel',[],...
    'ytick',[],'yticklabel',[],'XColor','white','YColor','white') 
xlabel('$r$ (arcsec)','Fontsize',f1,'Interpreter', 'Latex','Color','black')

subplot('position',[0.08 0.5 0.01 0.1])
set(gca,'FontSize',f1,...
    'TickLength', [0.02 0.02],'XMinorTick','off',...
    'YMinorTick','off','xtick',[],'xticklabel',[],...
    'ytick',[],'yticklabel',[],'XColor','white','YColor','white') 
ylabel('Velocity (km s$^{-1}$)','Fontsize',f1,'Interpreter', 'Latex','Color','black')

subplot('position',[0.16 0.13 0.83 0.86])

for i=1:length(v)-11
     top=v(i,1)+v_ci(i,1);
     bott=v(i,1)-v_ci(i,1);
     fill([r_ring(i,1), r_ring(i,2), r_ring(i,2), r_ring(i,1)],....
        [top top bott bott],[1,0.80,0.80],....
        'edgecolor',[1,0.8,0.8])
    
    if i==1
        hold on
    end
end

for i=1:length(v)-11
    plot([r_ring(i,1) r_ring(i,2)],(v(i,1)+v_ci(i,1)).*ones(2,1),'-','Color',[0.7 0 0],'LineWidth',1)
    plot([r_ring(i,1) r_ring(i,2)],(v(i,1)-v_ci(i,1)).*ones(2,1),'-','Color',[0.7 0 0],'LineWidth',1)
end

hold off

% do not change axis
axis([-5 170 120 230])

set(gca,'FontSize',f2,...
    'TickLength', [0.02 0.02],'XMinorTick','on',...
    'YMinorTick','on','xtick',xtick1,'xticklabel',xtick_label1,...
    'ytick',ytick1,'yticklabel',ytick_label1)


a='C:\Users\nit_n\Documents\MATLAB\h_alpha\n4321\figs\finished\vel_fit.eps';
set(figure(1), 'color', 'white')
print(figure(1),'-r600', '-depsc', a)

