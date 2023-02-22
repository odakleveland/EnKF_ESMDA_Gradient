%% EnKF/ES-MDA for multiple reflection in time with wavelet. 
clear all;close all;clc;

%% Input:
%Time
dt=0.002;                               %Time sample
t_max=dt*140;                           %max time
time=linspace(0,t_max,t_max/dt);        %time vector

theta=0:4:40;                           %takeoff angles

% nr=[60,100]+20;                        %Interface in time step
% vpt=[2300,2350,2400];vst=[1170,1200,1250];rhot=[2146,2192,2234];
% vpt=[2300,2500,2800];vst=[1170,1270,1440];rhot=[2146,2192,2234];

nr=[30,45,65,80,85,95,110,115]+20; %Interface in time step
vpt=[2300,2500,2150,2250,2400,2500,2300,2250,2400];
vst=[1170,1270,1070,1120,1170,1270,1170,1120,1170];
rhot=[2146,2192,2135,2110,2169,2192,2146,2110,2169];

%Background velocity and density
bvpt=vpt(1:end-1);bvst=vst(1:end-1);brhot=rhot(1:end-1);

%EnKF or ES-MDA? EnKF -> num=0, ES-MDA -> num=1
num=1;

[vp1D,vs1D,rho1D] = vel_den_vectors(time,nr,vpt,vst,rhot); %Create vectors for p,s and density plot.

figure('Name','Velocity model')
simple_model_plotting(time,vp1D,vs1D,rho1D)


%% Forward modeling
[G,d,a_alpha,a_beta,a_rho] = lin_zoeppritz(vpt,vst,rhot,theta,nr);

% a = 8e-03;b=9e-03; %snr ca 5
% a = 4e-03;b=7e-03; %snr ca 7
a = 1e-03;b=4.5e-03; %snr ca 15
% a = 6e-04;b=6e-04; %snr ca 70
noise = a + (b-a).*rand(88,1); % må endre ift størrelse på modell
d_noise = d + noise;

SNR_d = rms(d)./rms(d-d_noise)

d_noise = d;

%% L-curve
I = eye(size(G'*G));

n_alpha = length(I);
alpha_val = logspace(-2,0.01,n_alpha);

for ii = 1:n_alpha
    alpha = alpha_val(ii).^2*I;
    m_est(ii,:) = inv(G'*G+alpha)*G'*d_noise;
    residual_norm(ii) = norm(G*(m_est(ii,:)')-d);
    solution_norm(ii) = norm(m_est(ii,:));
end

figure('Name','L-curve')
plot(residual_norm.^2,solution_norm.^2,'-x','LineWidth',2); 
grid on
xlabel('Residual norm ||Gm_{est} - d||_2^2')
ylabel('Solution norm ||m_{est}||_2^2')
title('L-curve')

%% Inversion
r = residual_norm.^2;
num_alpha = find(r==min(r)); %from L-curve
[vp_inv,vs_inv,rho_inv,m_inv,m_est] = det_inversion_damped(G,d_noise,nr,vpt,vst,rhot,num_alpha);
%[vp_inv, vs_inv, rho_inv] = det_inversion(G,d_noise,nr,vpt,vst,rhot);

%% Ricker Wavelet
t = -t_max/2+dt:dt:t_max/2-dt;  
centr_frq = 25; 

w = ricker_wavelet_zero_phased(dt,t,centr_frq);

%% Create wavelet 
d_new = d';
d_conv = zeros(1,length(time)); %endret så det er mulig å plotte mot tid, riktig?
for kk = 1:length(d)
    d_conv(kk)=conv(d_new(:,kk),w,'same');
end

A = [a_alpha,a_beta,a_rho];
%A*m_est; %passer ikke

%%
% figure('Renderer', 'painters', 'Position', [50 40 800 500])
% hold on,title('No noise'),grid on
% for kk = 1:length(d)
%     %plot()
% end

%% Output

[vp1D_inv,vs1D_inv,rho1D_inv] = vel_den_vectors(time,nr,vp_inv,vs_inv,rho_inv);

% figure('Name','Inverted model')
% subplot(1,3,1)
% stairs(vp1D_inv,time,'r','Linewidth',1.5),title('P-wave'),grid on
% hold on
% stairs(vp1D,time,'k','Linewidth',1.5),title('P-wave'),grid on
% hold off
% xlabel('Velocity (m/s)'),ylabel('Time (s)')
% set(gca,'Ydir','reverse'),set(gca,'FontSize',10),set(gca,'Linewidth',2)
% axis([1000,4200,0,max(time)])
% subplot(1,3,2)
% stairs(vs1D_inv,time,'r','Linewidth',1.5),title('S-wave'),grid on
% hold on
% stairs(vs1D,time,'k','Linewidth',1.5),title('S-wave'),grid on
% hold off
% xlabel('Velocity (m/s)'),ylabel('Time (s)')
% set(gca,'Ydir','reverse'),set(gca,'FontSize',10),set(gca,'Linewidth',2)
% axis([500,2600,0,max(time)])
% subplot(1,3,3)
% stairs(rho1D_inv,time,'r','Linewidth',1.5),title('Density'),grid on
% hold on
% stairs(rho1D,time,'k','Linewidth',1.5),title('Density'),grid on
% legend('inverted','original ')
% hold off
% xlabel('Density (kg/m^3)'),ylabel('Time (s)')
% set(gca,'Ydir','reverse'),set(gca,'FontSize',10),set(gca,'Linewidth',2)
% axis([1800,2500,0,max(time)])

%% Covariance matrix
I=3;                                %Parameters per layer

%Reflection coefficients                         
Rpp=zeros((length(vp1D)),length(theta));     
for nn=1:length(vp1D)-1
    for kk=1:length(theta)
        
        a=vp1D(nn);aa1=vp1D(nn+1)-vp1D(nn);
        b=vs1D(nn);bb1=vs1D(nn+1)-vs1D(nn);
        rho=rho1D(nn);rhoo1=rho1D(nn+1)-rho1D(nn);
        
        %Non-linear zoeppritz to model the observed reflection coeff. Noise free
        Rpp(nn,kk)=linear_zoeppritz_contrast(theta(kk),a,aa1,b,bb1,rho,rhoo1); 
        %Rpp(nn,kk)=zoeppritz(theta(kk),a,aa1,b,bb1,rho,rhoo1); 
    end
end

obs1 = reflectivity_convolution(time,theta,Rpp,w);

%Priori model, assume gaussian
mu=zeros(1,I*length(nr));

%Constant prior mean
for ii=1:length(nr)
    mup(ii)=2350;mus(ii)=1270;mur(ii)=2150;
end

%Linear prior mean
% mup=linspace(vpt(1),max(vpt),length(nr));
% mus=linspace(vst(1),max(vst),length(nr));
% mur=linspace(rhot(1),max(rhot),length(nr));

%Gather the prior mean in one vector
for ii=1:length(nr)
    mu(1+I*(ii-1))=mup(ii);mu(2+I*(ii-1))=mus(ii);mu(3+I*(ii-1))=mur(ii);
end

Num_parameters=length(mu);      %Number of parameters in total
Pvari=500;Svari=300;Rvari=300;  %variance
Covariance_matrix=eye(Num_parameters,Num_parameters); 
for ii=1:length(nr)
    Covariance_matrix(1+I*(ii-1),1+I*(ii-1))=Pvari;
    Covariance_matrix(2+I*(ii-1),2+I*(ii-1))=Svari;
    Covariance_matrix(3+I*(ii-1),3+I*(ii-1))=Rvari;
    
    Covariance_matrix(2+I*(ii-1),1+I*(ii-1))=0; Covariance_matrix(1+I*(ii-1),2+I*(ii-1))=0;
    Covariance_matrix(3+I*(ii-1),1+I*(ii-1))=0;  Covariance_matrix(1+I*(ii-1),3+I*(ii-1))=0; 
    Covariance_matrix(3+I*(ii-1),2+I*(ii-1))=0;  Covariance_matrix(2+I*(ii-1),3+I*(ii-1))=0; 
end
Num_ensembles=1000;
[X]=ensemble_correlation(mu,Covariance_matrix,Num_parameters,Num_ensembles,1);   %Function to sample ensembles

%Add noise to noise free data
rng(41)                                         %Seed number, not needed
S_N_r=linspace(SNR_d,SNR_d,length(time));             %SNR, endret fra 15 til 100
%obs=zeros(1,length(time)*length(theta));
timepost=zeros(1,length(time)*length(theta));
for kk=1:length(theta)
    varR=rms(obs1((kk-1)*length(time)+1:kk*length(time)))./(S_N_r);     %STD
    %Observed data with noise
    obs((kk-1)*length(time)+1:kk*length(time))=...
        obs1((kk-1)*length(time)+1:kk*length(time))+randn(1,length(time)).*varR;
end

diagR=eye(1,length(obs1));
for kk=1:length(theta)
    diagR((kk-1)*length(time)+1:kk*length(time))=varR.^2;
end
R=diag(diagR);          %Error covariance matrix

plot_wavelets(obs,obs1,time,theta)

SNR=rms(obs1)./rms(obs1-obs)

%Change R - fix this later
diagRd=eye(1,length(vpt)-1);
for kk=1:length(theta)
    varR=rms(d((kk-1)*(length(vpt)-1)+1:kk*(length(vpt)-1)))./(SNR_d);     %STD
    diagRd((kk-1)*(length(vpt)-1)+1:kk*(length(vpt)-1))=varR.^2;
end
Rd=diag(diagRd);          %Error covariance matrix

%% Inversion part
tic

%ES-MDA
if num==1
    Na=5;
end
%EnKF
if num==0;
    Na=1;
end

a=[100/5,100/10,100/15,100/30,100/40]; %inflation coefficient
Xp=ESMDA_wavelet_linear(X,w,a,obs,R,Num_ensembles,Num_parameters,Na,bvpt,bvst,brhot,theta,time,nr);

time_toc=toc
%% Gradient-based method

I = 3; %number of parameters

Num_parameters=length(nr)*I;      %Number of parameters in total
cov_m=eye(Num_parameters); %covariance model
for ii=1:length(nr)
    cov_m(1+I*(ii-1),1+I*(ii-1))=Pvari;
    cov_m(2+I*(ii-1),2+I*(ii-1))=Svari;
    cov_m(3+I*(ii-1),3+I*(ii-1))=Rvari;
    
    cov_m(2+I*(ii-1),1+I*(ii-1))=0; cov_m(1+I*(ii-1),2+I*(ii-1))=0;
    cov_m(3+I*(ii-1),1+I*(ii-1))=0;  cov_m(1+I*(ii-1),3+I*(ii-1))=0; 
    cov_m(3+I*(ii-1),2+I*(ii-1))=0;  cov_m(2+I*(ii-1),3+I*(ii-1))=0; 
end

cov_d = eye(length(d))*max(max(Rd)); %covariance data

mean_posterior = m_est+cov_m*G'*inv(G*cov_m*G'+cov_d)*(d_noise-G*m_est); % eq 3.37 tarantola
m_posterior = inv(G'*inv(cov_d)*G+inv(cov_m));
conf_area = diag(m_posterior);conf_area_p = conf_area(1:3:end);conf_area_s = conf_area(2:3:end);conf_area_d = conf_area(3:3:end);

%% change from contrast to absolute velocities
mean_posterior_vel = zeros(length(nr),3);conf_area_vel_max = zeros(length(nr),3);conf_area_vel_min = zeros(length(nr),3);

for i = 1:length(nr)
    mean_posterior_vel(i,1:3) = [mean_posterior(3*(i-1)+1)*vpt(i)+vpt(i),mean_posterior(3*(i-1)+2)*vst(i)+vst(i),mean_posterior(3*i)*rhot(i)+rhot(i)];
end

%create vectors for velocities after inversion
vp_mean = [vpt(1), mean_posterior_vel(:,1)'];vs_mean = [vst(1), mean_posterior_vel(:,2)'];rho_mean = [rhot(1), mean_posterior_vel(:,3)'];
for i = 1:length(nr)
    conf_area_vel_max(i,1:3) = [conf_area_p(i)*vp_mean(i+1)+vp_mean(i+1),conf_area_s(i)*vs_mean(i+1)+vs_mean(i+1),conf_area_d(i)*rho_mean(i+1)+rho_mean(i+1)];
    conf_area_vel_min(i,1:3) = [-conf_area_p(i)*vp_mean(i+1)+vp_mean(i+1),-conf_area_s(i)*vs_mean(i+1)+vs_mean(i+1),-conf_area_d(i)*rho_mean(i+1)+rho_mean(i+1)];
end

vp_min = [vpt(1),conf_area_vel_min(:,1)'];vs_min = [vst(1),conf_area_vel_min(:,2)'];rho_min = [rhot(1),conf_area_vel_min(:,3)'];
vp_max = [vpt(1),conf_area_vel_max(:,1)'];vs_max= [vst(1),conf_area_vel_max(:,2)'];rho_max = [rhot(1),conf_area_vel_max(:,3)'];

%make new length for plotting
[vp1D_mean,vs1D_mean,rho1D_mean] = vel_den_vectors(time,nr,vp_mean,vs_mean,rho_mean);
[vp1D_min,vs1D_min,rho1D_min] = vel_den_vectors(time,nr,vp_min,vs_min,rho_min);
[vp1D_max,vs1D_max,rho1D_max] = vel_den_vectors(time,nr,vp_max,vs_max,rho_max);

%% Plotting fill:

figure(4)
plot_nl_wavelet_ESMDA_gradient_fill(X,Xp,vpt,vst,rhot,time,nr,I,1,vp1D_max,vp1D_min,vs1D_max,vs1D_min,rho1D_max,rho1D_min,vp1D,vs1D,rho1D)

%% Plotting mean:

figure(5)
plot_nl_wavelet_ESMDA_gradient(X,Xp,vpt,vst,rhot,time,nr,I,1,vp1D_mean,vs1D_mean,rho1D_mean)

% %% Plot gradient based method:
% 
% figure(6)
% plot_gradient(vp1D,vs1D,rho1D,vp1D_min,vp1D_max,vs1D_min,vs1D_max,rho1D_min,rho1D_max,time,vp1D_mean,vs1D_mean,rho1D_mean,nr)
% 
% %% Plot ES-MDA
% 
% figure(7)
% plot_nl_wavelet_ESMDA(X,Xp,vpt,vst,rhot,time,nr,I,1);
