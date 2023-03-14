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

%% Linear zoeppritz, not time dependent

delta_alpha = diff(vpt);

for i = 1:length(nr)
    for j = 1:length(theta)
        a_alpha(j,i) = 0.5*(1+(tand(theta(j)).^2));
    end  
end
a_alpha = reshape(a_alpha,[],1);

G = zeros(length(theta)*length(nr),length(nr));

for i = 1:length(nr)
    G(((i-1)*length(theta))+1:i*length(theta),((i-1))+1:i) = [a_alpha((i-1)*11+1:i*length(theta))];
end

%% Linear zoeppritz, time dependent

delta_alpha = diff(vp1D);delta_alpha(end+1)=0;

for i = 1:length(time)
    for  j = 1:length(theta)
        a_alpha(j,i) = 0.5*(1+(tand(theta(j)).^2));
        d(j,i) = a_alpha(j)*(delta_alpha(i)/vp1D(i));
    end
end

%% Convolution

obs_gradient = reflectivity_convolution(time,theta,d',w);

%% Adding noise to convolved data

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

%% Extract reflection coefficients from obs_gradient

obs_gradient_matrix = reshape(obs_gradient_noise,140,11);

d_obs = zeros(8,11);

for i = 1:length(nr)
    d_obs(i,:)=obs_gradient_matrix(nr(i),:);
end

d_obs = reshape(d_obs',88,1);


%% Sjekker om det blir perfekt resultat for d

d_obs_test = nonzeros(d);

%% Covariance d

diagRd=eye(1,length(vpt)-1);
for kk=1:length(theta)
    varR=rms(d_obs((kk-1)*(length(vpt)-1)+1:kk*(length(vpt)-1)))./(max(S_N_r));     %STD
    diagRd((kk-1)*(length(vpt)-1)+1:kk*(length(vpt)-1))=varR.^2;
end
cov_d=diag(diagRd);          %Error covariance matrix

%% Gradient based method

I = 1; %number of parameters
Pvari=500;
Num_parameters=length(nr)*I;      %Number of parameters in total *I
cov_m=eye(Num_parameters)*Pvari; %140x140
%cov_d = eye(length(theta)*length(nr))*max(max(R)); %kan være en feil her, se hovedscript

m_prior = zeros(1,length(nr))*2350;

obs_gradient_matrix = reshape(obs_gradient,140,11);

%Endre: bruk riktig verdi, ikke vpt!!!
part1 = cov_m*G';
part2 = inv(G*cov_m*G'+cov_d);
part3  = d_obs-(G*m_prior');

mean_posterior = m_prior'+part1*part2*part3;

%% Gradient based method with convolution of G

%does not work

% G = a_alpha(:,1);
% 
% G_conv=conv(G,w,'same');
% 
% % part1 = cov_m*G_conv;
% % part2 = inv(G_conv'*cov_m*G_conv+cov_d);
% % part3  = obs_gradient_matrix-(G_conv*m_est)';
% part1 = G_conv;
% part2 = inv(G_conv'*G_conv);
% part3  = obs_gradient_matrix-(G_conv*m_prior)';
% 
% test = m_prior+part1*part2*part3';

%% Change from contrast to velocity

m_post_inv = zeros(length(nr),1);

for i = 1:length(nr)
    m_post_inv(i) = mean_posterior(i)*vpt(i)+vpt(i);
end

m_post_inv = [vpt(1); m_post_inv];

%% Confidence area

%not working, but works in SC_ESMDA_G for p-wave, s-wave and density

m_posterior = inv(G'*inv(cov_d)*G+inv(cov_m));
conf_area = diag(m_posterior);

conf_area_max= zeros(1,length(nr));conf_area_min= zeros(1,length(nr));

for i = 1:length(nr)
    conf_area_max(i) = conf_area(i)*m_post_inv(i+1)+m_post_inv(i+1);
    conf_area_min(i) = -conf_area(i)*m_post_inv(i+1)+m_post_inv(i+1);
end
conf_area_max = [vpt(1);conf_area_max'];conf_area_min = [vpt(1);conf_area_min'];

%% change length of vectors to length(time)

nnr=[0,nr,length(time)];
for gg=1:length(nnr)-1
    vp1D_mean(1+nnr(gg):nnr(gg+1))=m_post_inv(gg);
    vp1D_max(1+nnr(gg):nnr(gg+1))=conf_area_max(gg);
    vp1D_min(1+nnr(gg):nnr(gg+1))=conf_area_min(gg);
end

%% Plot velocities
%set axis values
x_min_vp=min(vp1D_min)-200*1.96;x_max_vp=max(vp1D_max)+300*1.96;
y_min=time(nr(1)+1);y_max=max(time);

figure(10)
%subplot(1,3,1)
hold all
p1 = patch([vp1D_max fliplr(vp1D_min)], [time fliplr(time)],[0.2,0.2,0.2]); set(p1, 'facealpha',.5)
%p1 = patch([vp1D_max fliplr(vp1D_min)], [time fliplr(time)],'r'); set(p1, 'facealpha',.5)
stairs(vp1D_mean,time,'r','Linewidth',3.5)
stairs(vp1D,time,'k','Linewidth',3.5),title('P-wave velocity'),grid on
hold off
xlabel('Velocity (m/s)'),ylabel('TWT (s)')
set(gca,'Ydir','reverse'),set(gca,'FontSize',10),set(gca,'Linewidth',2)
legend('95% posterior confidence area','Posterior mean','True model')
axis([x_min_vp,x_max_vp,y_min,y_max]),set(gca,'FontSize',35)



