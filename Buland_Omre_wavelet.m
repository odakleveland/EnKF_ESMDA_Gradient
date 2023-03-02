clear;close all;clc

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

new_m_est = zeros(9,3);
for i = 0:2
    new_m_est(2:9,i+1) = m_est(i*8+1:8*i+8);
end
[m_est_alpha,m_est_beta,m_est_rho] = vel_den_vectors(time,nr,new_m_est(:,1),new_m_est(:,2),new_m_est(:,3));

m_est = [m_est_alpha',m_est_beta',m_est_rho'];
%% Finding G and d that depends on time

[G_new,d_new,a_alpha,a_beta,a_rho] = lin_zoeppritz_new(vp1D,vs1D,rho1D,theta,nr,time);

% [vp_inv_new,vs_inv,rho_inv,m_inv,m_est] = det_inversion_damped(G,d,nr,vpt,vst,rhot,1);

%%

% for i = 1:length(time)
%     for j = 1:length(theta)
%         a_alpha(i,j) = 0.5*(1+(tand(theta(j)).^2)); %må få denne til å avhenge av tid
%     end  
% end

%%
m_est_alpha = [m_est(:,1)];

c = a_alpha.*m_est_alpha;

%%

d_obs = zeros(length(time),length(theta));
for kk = 1:length(theta)
    con = conv(c(:,kk),w,'same');
    d_obs(:,kk) = con;
end

%% plot wavelet
normm=norm(d_obs(1:length(time)));
figure('Renderer', 'painters', 'Position', [50 40 800 500])
for kk = 1:length(theta)
    plot(-d_obs(:,kk)/normm*7*4+((kk-1)*4),time)
    hold on
end
stairs(vp1D*10^(-03),time,'k','Linewidth',1.5),title('P-wave'),grid on
hold off

%% test

% delta_alpha = diff(vp1D);delta_beta = diff(vs1D);delta_rho = diff(rho1D);
% delta_alpha(end+1)=0;delta_beta(end+1)=0;delta_rho(end+1)=0;
% 
% G = zeros(length(time),3*length(nr));
% d = zeros(length(time),length(theta));
% 
% for i = 1:length(time)
%     for j = 1:length(theta)
%         a_alpha(i,j) = 0.5*(1+(tand(theta(j)).^2));
%         a_beta(i,j) = -4*(vs1D(i).^2)/(vp1D(i).^2)*sind(theta(j)).^2;
%         a_rho(i,j) = 0.5*(1-4*(vs1D(i).^2)/(vp1D(i).^2)*sind(theta(j)).^2);
%         d(i,j) = a_alpha(i,j)*(delta_alpha(i)/vp1D(i))+a_beta(i,j)*(delta_beta(i)/vs1D(i))+a_rho(i,j)*(delta_rho(i)/rho1D(i));
%     end  
% end
% 
% G = [a_alpha,a_beta,a_rho];
% 
% %d = [d;d;d];
% 
% size(d)
% size(G)
% 
% num_alpha=1;
% I = eye(size(G'*G));
% 
% n_alpha = length(I);
% alpha = logspace(-2,0.01,n_alpha);
% 
% m_est = inv(G'*G+alpha(num_alpha).^2*I)*(G'*d);
% 
% size(a_alpha)

%% Buland + Omre eq. 22 (inverse modeling)
% want to find m from d = S*Am
