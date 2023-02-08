function [] = simple_model_plotting(time,vp1D,vs1D,rho1D)

%plot hastighetsmodell
subplot(1,3,1)
stairs(vp1D,time,'k','Linewidth',1.5),title('P-wave'),grid on
xlabel('Velocity (m/s)'),ylabel('Time (s)')
set(gca,'Ydir','reverse'),set(gca,'FontSize',10),set(gca,'Linewidth',2)
axis([1000,4200,0,max(time)])
subplot(1,3,2)
stairs(vs1D,time,'k','Linewidth',1.5),title('S-wave'),grid on
xlabel('Velocity (m/s)'),ylabel('Time (s)')
set(gca,'Ydir','reverse'),set(gca,'FontSize',10),set(gca,'Linewidth',2)
axis([500,2600,0,max(time)])
subplot(1,3,3)
stairs(rho1D,time,'k','Linewidth',1.5),title('Density'),grid on
xlabel('Density (kg/m^3)'),ylabel('Time (s)')
set(gca,'Ydir','reverse'),set(gca,'FontSize',10),set(gca,'Linewidth',2)
axis([1800,2500,0,max(time)])

end