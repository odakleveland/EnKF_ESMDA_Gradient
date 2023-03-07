function plot_gradient(vp1D,vs1D,rho1D,vp1D_min,vp1D_max,vs1D_min,vs1D_max,rho1D_min,rho1D_max,time,vp1D_mean,vs1D_mean,rho1D_mean,nr)
%Set axis values:
x_min_vp=min(vp1D_min)-200*1.96;x_max_vp=max(vp1D_max)+300*1.96;
x_min_vs=min(vs1D_min)-200*1.96;x_max_vs=max(vs1D_max)+200*1.96;
x_min_rho=min(rho1D_min)-200*1.96;x_max_rho=max(rho1D_max)+200*1.96;
y_min=time(nr(1)+1);y_max=max(time);

%figure('Name','Gradient based methods')
subplot(1,3,1)
hold all
p1 = patch([vp1D_max fliplr(vp1D_min)], [time fliplr(time)],[0.2,0.2,0.2]); set(p1, 'facealpha',.5)
%p1 = patch([vp1D_max fliplr(vp1D_min)], [time fliplr(time)],'r'); set(p1, 'facealpha',.5)
stairs(vp1D_mean,time,'r','Linewidth',3.5)
stairs(vp1D,time,'k','Linewidth',3.5),title('P-wave velocity'),grid on
hold off
xlabel('Velocity (m/s)'),ylabel('TWT (s)')
set(gca,'Ydir','reverse'),set(gca,'FontSize',10),set(gca,'Linewidth',2)
axis([x_min_vp,x_max_vp,y_min,y_max]),set(gca,'FontSize',35)

subplot(1,3,2)
hold all
p2 = patch([vs1D_max fliplr(vs1D_min)], [time fliplr(time)], [0.2,0.2,0.2]); set(p2, 'facealpha',.5)
stairs(vs1D_mean,time,'r','Linewidth',3.5)
stairs(vs1D,time,'k','Linewidth',3.5),title('S-wave velocity'),grid on
hold off
xlabel('Velocity (m/s)'),ylabel('TWT (s)')
set(gca,'Ydir','reverse'),set(gca,'FontSize',10),set(gca,'Linewidth',2)
axis([x_min_vs,x_max_vs,y_min,y_max]),set(gca,'FontSize',35)

subplot(1,3,3)
hold all
p3 = patch([rho1D_max fliplr(rho1D_min)], [time fliplr(time)],[0.2,0.2,0.2]); set(p3, 'facealpha',.5)
stairs(rho1D_mean,time,'r','Linewidth',3.5)
stairs(rho1D,time,'k','Linewidth',3.5),title('Density'),grid on
legend('95% posterior confidence area','Posterior mean','True model')
hold off
xlabel('Density (kg/m^3)'),ylabel('TWT (s)')
set(gca,'Ydir','reverse'),set(gca,'FontSize',10),set(gca,'Linewidth',2)
axis([x_min_rho,x_max_rho,y_min,y_max]),set(gca,'FontSize',35)
end