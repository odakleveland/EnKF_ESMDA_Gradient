function simple_model_inverse_plotting(vp1D_inv,vp1D,vs1D_inv,vs1D,rho1D_inv,rho1D,time)

figure('Name','Inverted model')
subplot(1,3,1)
stairs(vp1D_inv,time,'r','Linewidth',1.5),title('P-wave'),grid on
hold on
stairs(vp1D,time,'k','Linewidth',1.5),title('P-wave'),grid on
hold off
xlabel('Velocity (m/s)'),ylabel('Time (s)')
set(gca,'Ydir','reverse'),set(gca,'FontSize',10),set(gca,'Linewidth',2)
axis([1000,4200,0,max(time)])
subplot(1,3,2)
stairs(vs1D_inv,time,'r','Linewidth',1.5),title('S-wave'),grid on
hold on
stairs(vs1D,time,'k','Linewidth',1.5),title('S-wave'),grid on
hold off
xlabel('Velocity (m/s)'),ylabel('Time (s)')
set(gca,'Ydir','reverse'),set(gca,'FontSize',10),set(gca,'Linewidth',2)
axis([500,2600,0,max(time)])
subplot(1,3,3)
stairs(rho1D_inv,time,'r','Linewidth',1.5),title('Density'),grid on
hold on
stairs(rho1D,time,'k','Linewidth',1.5),title('Density'),grid on
legend('inverted','original ')
hold off
xlabel('Density (kg/m^3)'),ylabel('Time (s)')
set(gca,'Ydir','reverse'),set(gca,'FontSize',10),set(gca,'Linewidth',2)
axis([1800,2500,0,max(time)])

end