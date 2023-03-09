clear;close all;clc;

%% Velocity model

dt=0.002;                               %Time sample
t_max=dt*140;                           %max time
time=linspace(0,t_max,t_max/dt);        %time vector

theta=0:4:40;                           %takeoff angles

nr=[30,45,65,80,85,95,110,115]+20; %Interface in time step
vpt=[2300,2500,2150,2250,2400,2500,2300,2250,2400]; %p-wave velocity

nnr=[0,nr,length(time)];
for gg=1:length(nnr)-1
    vp1D(1+nnr(gg):nnr(gg+1))=vpt(gg);
end

%% Ricker wavelet

t = -t_max/2+dt:dt:t_max/2-dt;  
centr_frq = 25; 

w = ricker_wavelet_zero_phased(dt,t,centr_frq);

%% Linear zoeppritz, time dependent

delta_alpha = diff(vp1D);delta_alpha(end+1)=0;

G = zeros(length(time),length(nr));

for i = 1:length(time)
    for  j = 1:length(theta)
        a_alpha(j,i) = 0.5*(1+(tand(theta(j)).^2));
        d(j,i) = a_alpha(j)*(delta_alpha(i)/vp1D(i));
    end
end

G = a_alpha(1:11,1);

%% Deterministic inversion

I = eye(size(G'*G));

n_alpha = length(I);
alpha = logspace(-2,0.01,n_alpha);

num_alpha=1;

m_est = inv(G'*G+alpha(num_alpha).^2*I)*G'*d;

%% Sjekk at m_est stemmer med hastigheter

m_est_test = nonzeros(m_est);

m_inv = zeros(length(nr),1);

for i = 1:length(nr)
    m_inv(i) = m_est_test(i)*vpt(i)+vpt(i);
end

%% Convolution

obs_gradient = reflectivity_convolution(time,theta,d',w);

%% Adding noise

S_N_r=linspace(7,7,length(time));             %SNR, endret fra 15 til 100
timepost=zeros(1,length(time)*length(theta));
for kk=1:length(theta)
    varR=rms(obs_gradient((kk-1)*length(time)+1:kk*length(time)))./(S_N_r);     %STD
    %Observed data with noise
    obs_gradient_noise((kk-1)*length(time)+1:kk*length(time))=...
        obs_gradient((kk-1)*length(time)+1:kk*length(time))+randn(1,length(time)).*varR;
end

SNR=rms(obs_gradient)./rms(obs_gradient-obs_gradient_noise)

diagR=eye(1,length(obs_gradient));
for kk=1:length(theta)
    varR=rms(obs_gradient((kk-1)*length(time)+1:kk*length(time)))./S_N_r;   %STD
    diagR((kk-1)*length(time)+1:kk*length(time))=varR.^2;
end
R=diag(diagR);

%% Plotting waveforms
plot_wavelets(obs_gradient_noise,obs_gradient,time,theta)

%% Gradient based method
%g-matrise trenger ikke avhenge av tid
%se på buland omre eq 30+31, sammenlign med tarantola

I = 1; %number of parameters
Pvari=500;Svari=300;Rvari=300;
Num_parameters=length(theta)*I;      %Number of parameters in total *I
cov_m=eye(Num_parameters)*Pvari;
cov_d = eye(length(1))*max(max(R));

obs_gradient_matrix = reshape(obs_gradient,140,11);

partt1 = cov_m*G;
partt2 = inv(G'*cov_m*G+cov_d);
partt3  = obs_gradient_matrix-(G*m_est)';
%partt3 = zeros(140,1); %test for å se om partt1 og partt2 er riktig

mean_posterior = m_est+partt1*partt2.*partt3';
mean_posterior = mean_posterior';
mean_posterior = mean_posterior(:,1);

%% Change from contrast to velocity

m_post_inv = zeros(length(time),1);

for i = 2:length(nr)
    m_post_inv(1:nr(1)) = vpt(1);
    m_post_inv(nr(i-1)+1:nr(i)) = mean_posterior(nr(i-1)+1)*vpt(i-1)+vpt(i);
    m_post_inv(nr(end)+1:end) = mean_posterior(nr(end-1)+1)*vpt(end-1)+vpt(end); %sjekk denne
end

%% Confidence area

m_posterior = inv(G'*inv(cov_d)*G+inv(cov_m));

%% Plot velocities

figure(2)
stairs(m_post_inv,time)
hold on
stairs(vp1D,time)
legend('inverted model','true model')


