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
m = [vpt',vst',rhot'];

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

%d_noise = d;

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

%% Output

[vp1D_inv,vs1D_inv,rho1D_inv] = vel_den_vectors(time,nr,vp_inv,vs_inv,rho_inv);

simple_model_inverse_plotting(vp1D_inv,vp1D,vs1D_inv,vs1D,rho1D_inv,rho1D,time)

%% Covariance matrix
I=3;                                %Parameters per layer

%Reflection coefficients                         
Rpp = reflection_coefficients(vp1D,vs1D,rho1D,theta);

%Convolution
obs1 = reflectivity_convolution(time,theta,Rpp,w);

%Priori model, assume gaussian
mu=zeros(1,I*length(nr));

%Constant prior mean
for ii=1:length(nr)
    mup(ii)=2350;mus(ii)=1270;mur(ii)=2150;
end

%Gather the prior mean in one vector
for ii=1:length(nr)
    mu(1+I*(ii-1))=mup(ii);mu(2+I*(ii-1))=mus(ii);mu(3+I*(ii-1))=mur(ii);
end

Num_parameters=length(mu);      %Number of parameters in total
Pvari=500;Svari=300;Rvari=300;  %variance

cov_m = covariance_matrix_m(I,Num_parameters,nr,Pvari,Svari,Rvari);

Num_ensembles=1000;
[X]=ensemble_correlation(mu,cov_m,Num_parameters,Num_ensembles,1);   %Function to sample ensembles

%Add noise to noise free data
rng(41)                                         %Seed number, not needed
S_N_r=linspace(SNR_d,SNR_d,length(time));             %SNR, endret fra 15 til 100
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

cov_d = eye(length(d))*max(max(Rd)); %covariance data
%cov_m = covariance_matrix_m(I,Num_parameters,nr,Pvari,Svari,Rvari)

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

%% Buland and Omre method (only P-waves)

%changing the size om m_est:
m_est_alpha =zeros(9,1);m_est_beta =zeros(9,1);m_est_rho =zeros(9,1);
for i = 1:length(nr)
    m_est_alpha(i+1) = m_est(3*(i-1)+1);
    m_est_beta(i+1) = m_est(3*(i-1)+2);
    m_est_rho(i+1) = m_est(3*i);
end

nnr=[0,nr,length(time)];
for gg=1:length(nnr)-1
    m_est_1D_alpha(1+nnr(gg):nnr(gg+1))=m_est_alpha(gg);
    m_est_1D_beta(1+nnr(gg):nnr(gg+1))=m_est_beta(gg);
    m_est_1D_rho(1+nnr(gg):nnr(gg+1))=m_est_rho(gg);
end

m_est_1D = [m_est_1D_alpha',m_est_1D_beta',m_est_1D_rho'];

% Finding G and d that depends on time

[G_new,d_new,a_alpha,a_beta,a_rho] = lin_zoeppritz_new(vp1D,vs1D,rho1D,theta,nr,time);

%Finding new m_inv

m_inv_alpha_1D=zeros(140,1);

m_inv_alpha1D = zeros(length(time),1);
for i = 2:length(nr)
    m_inv_alpha_1D(1:nr(1)) = vpt(1);
    m_inv_alpha_1D(nr(i-1)+1:nr(i)) = m_est_1D_alpha(nr(i))*vpt(i-1)+vpt(i-1);
    m_inv_alpha_1D(nr(end):end) = m_est_1D_alpha(nr(end)+1)*vpt(end-1)+vpt(end-1);
end

%Convolution
obs_d_gradient = reflectivity_convolution(time,theta,d_new,w);

% Adding noise
for kk=1:length(theta)
    varR=rms(obs_d_gradient((kk-1)*length(time)+1:kk*length(time)))./(S_N_r);     %STD
    %Observed data with noise
    obs_d_gradient_noise((kk-1)*length(time)+1:kk*length(time))=...
        obs_d_gradient((kk-1)*length(time)+1:kk*length(time))+randn(1,length(time)).*varR;
end

% Plotting waveforms
plot_wavelets(obs_d_gradient_noise,obs_d_gradient,time,theta)
SNR_gradient=rms(obs_d_gradient)./rms(obs_d_gradient-obs_d_gradient_noise)

%% Gradient based method, only p-waves

I_new = 1; %number of parameters
Pvari=500;Svari=300;Rvari=300;
Num_parameters=length(theta)*I_new;      %Number of parameters in total *I
cov_m_new=eye(Num_parameters)*Pvari; %covariance model
cov_d_new = eye(length(d_new))*max(max(Rd));

G_new_alpha=G_new(:,1:11);

obs_d_gradient_matrix = reshape(obs_d_gradient_noise,140,11);

partt1 = cov_m_new*G_new_alpha';
partt2 = inv(G_new_alpha*cov_m_new*G_new_alpha'+cov_d_new);
partt3  = obs_d_gradient_matrix-G_new_alpha.*m_est_1D(:,1);

mean_posterior = m_est_1D_alpha+partt1*partt2.*partt3';
mean_posterior=mean_posterior';

%change from contrast to velocities

mean_posterior_vel_test = zeros(length(time),1);
    
for i = 2:length(nr)
    mean_posterior_vel_test(1:nr(1)) = vpt(1);
    mean_posterior_vel_test(nr(i-1)+1:nr(i)) = mean_posterior(nr(i))*vpt(i-1)+vpt(i-1);
    mean_posterior_vel_test(nr(end)+1:end) = mean_posterior(nr(end)+1)*vpt(end-1)+vpt(end-1);
end

%% Plotting fill:

figure(4)
plot_nl_wavelet_ESMDA_gradient_fill(X,Xp,vpt,vst,rhot,time,nr,I,1,vp1D_max,vp1D_min,vs1D_max,vs1D_min,rho1D_max,rho1D_min,vp1D,vs1D,rho1D)

%% Plotting mean:

figure(5)
plot_nl_wavelet_ESMDA_gradient(X,Xp,vpt,vst,rhot,time,nr,I,1,vp1D_mean,vs1D_mean,rho1D_mean)

%% Plot gradient based method:

figure(6)
plot_gradient(vp1D,vs1D,rho1D,vp1D_min,vp1D_max,vs1D_min,vs1D_max,rho1D_min,rho1D_max,time,mean_posterior_vel_test,vs1D_mean,rho1D_mean,nr)

%% Plot ES-MDA

figure(7)
plot_nl_wavelet_ESMDA(X,Xp,vpt,vst,rhot,time,nr,I,1);

%% Test


[G_new,d_new,a_alpha,a_beta,a_rho] = lin_zoeppritz_new(vp1D,vs1D,rho1D,theta,nr,time);
