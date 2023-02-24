clear all;close all;clc

%% Velocity model

dt=0.002;                               %Time sample
t_max=dt*140;                           %max time
time=linspace(0,t_max,t_max/dt);        %time vector

theta=0:4:40;                           %takeoff angles

% nr=[40]+20;                        %Interface in time step
% vpt=[2300,2450];vst=vpt*0.5;rhot=[2146,2192];
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

A = [a_alpha,a_beta,a_rho];
S = w;
c = G*m_est;

d_obs = zeros(1,length(d));
for kk = 1:length(d)
    d_obs(:,kk)=conv(c(kk,:),S,'same');
end
max(max(d-d_obs'))

%% Buland + Omre eq. 22 (inverse modeling)
% want to find m from d = S*Am


