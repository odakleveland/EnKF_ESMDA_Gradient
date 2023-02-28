clear all;close all;clc

%% Velocity model

dt=0.002;                               %Time sample
t_max=dt*140;                           %max time
time=linspace(0,t_max,t_max/dt);        %time vector

theta=0:4:40;                           %takeoff angles

% nr=[30,45]+20;                        %Interface in time step
% vpt=[2300,2450,2150];vst=vpt*0.5;rhot=[2146,2192,2135];
nr=[30,45,65,80,85,95,110,115]+20; %Interface in time step
vpt=[2300,2500,2150,2250,2400,2500,2300,2250,2400];
vst=[1170,1270,1070,1120,1170,1270,1170,1120,1170];
rhot=[2146,2192,2135,2110,2169,2192,2146,2110,2169];

[vp1D,vs1D,rho1D] = vel_den_vectors(time,nr,vpt,vst,rhot); %Create vectors for p,s and density plot.

figure('Name','Velocity model')
simple_model_plotting(time,vp1D,vs1D,rho1D)

%% Ricker wavelet

t = -t_max/2+dt:dt:t_max/2-dt;  
centr_frq = 25; 

w = ricker_wavelet_zero_phased(dt,t,centr_frq);

%% Linearized Zoeppritz (Aki & Richards)

[G,d,a_alpha,a_beta,a_rho] = lin_zoeppritz(vpt,vst,rhot,theta,nr);

%% Deterministic inversion to find m_est

[vp_inv,vs_inv,rho_inv,m_inv,m_est] = det_inversion_damped(G,d,nr,vpt,vst,rhot,1);


%% Buland + Omre eq. 22 (forward modeling)


S = w; %kan bli matrise (en wavelet per vinkel)
c = G*m_est; %skal avhenge av vinkel og tid

%test for 1 laggrense først
%endre på laggrensen 
% 
% d_obs=zeros(1,((length(time)-1)*length(theta)));
% for ii = 1:length(nr)
%     d_obs((kk-1)*length(time)+1:kk*length(time))=conv(Rpp(:,ii),S,'same');
% end

% d_obs = reflectivity_convolution(time,theta,Rpp,w);

% figure('Renderer', 'painters', 'Position', [50 40 800 500])
% subplot(1,3,1)
% normm=norm(d_obs(1:length(time)));
% hold on,title('With noise'),grid on
% for kk=1:length(theta)
%     plot(d_obs((kk-1)*length(time)+1:kk*length(time))/normm*7*4+((kk-1)*4),time(1:end),'k','Linewidth',2)
% end
% set(gca,'Ydir','reverse'),set(gca,'FontSize',15),set(gca,'Linewidth',2)
% xlabel('Incidence angle'),ylabel('Time (s)')
% hold off
% axis([-10,50,0.1,max(time)])

%% Buland Omre eq. 22 (forward modeling, p-waves)

m_est_alpha = [m_est(1:8);m_est(1:8);m_est(1:8);m_est(1:8);m_est(1:8);m_est(1:8);m_est(1:8);m_est(1:8);m_est(1:8);m_est(1:8);m_est(1:8)];

c = a_alpha*m_est_alpha';

d_obs = zeros(88,length(theta));
for kk = 1:length(theta)
    con = conv(c(:,kk),S,'same');
    d_obs(:,kk) = con;
end

%plot time(51:end-2) (første del av time er for første lag)
%% plot wavelet
figure('Renderer', 'painters', 'Position', [50 40 800 500])
for kk = 1:length(theta)
    plot(d_obs(:,kk),time(51:end-2))
    hold on
end
hold off
