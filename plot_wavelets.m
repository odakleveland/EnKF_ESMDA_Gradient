function plot_wavelets(obs,obs1,time,theta)

figure('Renderer', 'painters', 'Position', [50 40 800 500])
subplot(1,3,1)
normm=norm(obs(1:length(time)));
hold on,title('With noise'),grid on
for kk=1:length(theta)
    plot(obs((kk-1)*length(time)+1:kk*length(time))/normm*7*4+((kk-1)*4),time(1:end),'k','Linewidth',2)
end
set(gca,'Ydir','reverse'),set(gca,'FontSize',15),set(gca,'Linewidth',2)
xlabel('Incidence angle'),ylabel('Time (s)')
hold off
axis([-10,50,0.1,max(time)])

subplot(1,3,2)
normm=norm(obs(1:length(time)));
hold on,title('Without noise'),grid on
for kk=1:length(theta)
    plot(obs1((kk-1)*length(time)+1:kk*length(time))/normm*7*4+((kk-1)*4),time(1:end),'k','Linewidth',2)
end
set(gca,'Ydir','reverse'),set(gca,'FontSize',15),set(gca,'Linewidth',2)
xlabel('Incidence angle'),ylabel('Time (s)')
hold off
axis([-10,50,0.1,max(time)])

subplot(1,3,3)
normm=norm(obs(1:length(time)));
hold on,title('Difference'),grid on
for kk=1:length(theta)
    diffobs=obs((kk-1)*length(time)+1:kk*length(time))-obs1((kk-1)*length(time)+1:kk*length(time));
    plot(diffobs/normm*7*4+((kk-1)*4),time(1:end),'k','Linewidth',2)
end
set(gca,'Ydir','reverse'),set(gca,'FontSize',15),set(gca,'Linewidth',2)
xlabel('Incidence angle'),ylabel('Time (s)')
hold off
axis([-10,50,0.1,max(time)])

end