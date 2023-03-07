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

%plot(w,t)

%% Linearized Zoeppritz (Aki & Richards)

[G,d,a_alpha,a_beta,a_rho] = lin_zoeppritz(vpt,vst,rhot,theta,nr);

%% Deterministic inversion to find m_est

[vp_inv,vs_inv,rho_inv,m_inv,m_est] = det_inversion_damped(G,d,nr,vpt,vst,rhot,1);

[vp1D_inv,vs1D_inv,rho1D_inv] = vel_den_vectors(time,nr,vp_inv,vs_inv,rho_inv);

m1D_inv = [vp1D_inv',vs1D_inv',rho1D_inv'];


%% endre størrelse av m_est

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


%% Finding G, d and m that depends on time

[G_new,d_new,a_alpha,a_beta,a_rho] = lin_zoeppritz_new(vp1D,vs1D,rho1D,theta,nr,time);

I = eye(size(G_new'*G_new));

n_alpha = length(I);
alpha = logspace(-2,0.01,n_alpha);

num_alpha=1;

m_est_new = inv(G_new'*G_new+alpha(num_alpha).^2*I)*(G_new'*d_new);

%% Find new m_inv

nnr=[0,nr,length(time)];
for gg=1:length(nnr)-1
    m_est_1D_alpha(1+nnr(gg):nnr(gg+1))=m_est_alpha(gg);
end

m_inv_alpha_1D=zeros(140,1);

m_inv_alpha1D = zeros(length(time),1);
for i = 2:length(nr)
    m_inv_alpha_1D(1:nr(1)) = vpt(1);
    m_inv_alpha_1D(nr(i-1)+1:nr(i)) = m_est_1D_alpha(nr(i))*vpt(i-1)+vpt(i-1);
    m_inv_alpha_1D(nr(end):end) = m_est_1D_alpha(nr(end)+1)*vpt(end-1)+vpt(end-1);
end

%% Convolution

obs_d_gradient = reflectivity_convolution(time,theta,d_new,w);

%% Adding noise

S_N_r=linspace(7,7,length(time));             %SNR, endret fra 15 til 100
%obs=zeros(1,length(time)*length(theta));
timepost=zeros(1,length(time)*length(theta));
for kk=1:length(theta)
    varR=rms(obs_d_gradient((kk-1)*length(time)+1:kk*length(time)))./(S_N_r);     %STD
    %Observed data with noise
    obs_d_gradient_noise((kk-1)*length(time)+1:kk*length(time))=...
        obs_d_gradient((kk-1)*length(time)+1:kk*length(time))+randn(1,length(time)).*varR;
end

%% Plotting waveforms
plot_wavelets(obs_d_gradient_noise,obs_d_gradient,time,theta)

%% Gradient-based method
%må finne ut hvordan dette blir med ny G og d + regne ut convariance som
%stemmer (cov_d+cov_m)
%må finne ny m_est
%refner med: d_noise=obs_d_gradient_noise

%covariance model:
I = 3; %number of parameters
Pvari=500;Svari=300;Rvari=300;
Num_parameters=length(theta)*I;      %Number of parameters in total *I
cov_m=eye(Num_parameters); %covariance model
for ii=1:length(nr)
    cov_m(1+I*(ii-1),1+I*(ii-1))=Pvari;
    cov_m(2+I*(ii-1),2+I*(ii-1))=Svari;
    cov_m(3+I*(ii-1),3+I*(ii-1))=Rvari;
    
    cov_m(2+I*(ii-1),1+I*(ii-1))=0; cov_m(1+I*(ii-1),2+I*(ii-1))=0;
    cov_m(3+I*(ii-1),1+I*(ii-1))=0;  cov_m(1+I*(ii-1),3+I*(ii-1))=0; 
    cov_m(3+I*(ii-1),2+I*(ii-1))=0;  cov_m(2+I*(ii-1),3+I*(ii-1))=0; 
end

%covariance data:
SNR=rms(obs_d_gradient)./rms(obs_d_gradient-obs_d_gradient_noise)

%Change R - fix this later
diagRd=eye(1,length(vpt)-1);
for kk=1:length(theta)
    varR=rms(d((kk-1)*(length(vpt)-1)+1:kk*(length(vpt)-1)))./(SNR);     %STD
    diagRd((kk-1)*(length(vpt)-1)+1:kk*(length(vpt)-1))=varR.^2;
end
Rd=diag(diagRd); 
cov_d = eye(length(d_new))*max(max(Rd));

%mean_posterior = m_est+cov_m*G'*inv(G*cov_m*G'+cov_d)*(d_noise-G*m_est); % eq 3.37 tarantola
%p-waves:
% part1=cov_m*G_new(:,1:11)';
% part2=inv(G_new(:,1:11)*cov_m*G_new(:,1:11)'+cov_d);
% part3=(d_new-G_new(:,1:11).*m1D_inv(:,1));
% 
% mean_posterior = m1D_inv(:,1)+part1*part2*part3;

% mean_posterior = m1D_inv(:,1)+cov_m*G_new(:,1:11)'*inv(G_new(:,1:11)*cov_m*G_new(:,1:11)'+cov_d)*(d_new-G_new(:,1:11).*m1D_inv(:,1)); % eq 3.37 tarantola
m_posterior = inv(G_new'*inv(cov_d)*G_new+inv(cov_m));


mean_posterior = m_est_new+cov_m*G_new'*inv(G_new*cov_m*G_new'+cov_d)*(d_new-G_new*m_est_new); % eq 3.37 tarantola
%% Buland + Omre eq. 22 (inverse modeling)
% want to find m from d = S*Am

%% Gradient based method, p-waves

I = 1; %number of parameters
Pvari=500;Svari=300;Rvari=300;
Num_parameters=length(theta)*I;      %Number of parameters in total *I
Num_parameters = 11;
cov_m=eye(Num_parameters)*Pvari; %covariance model
% for ii=1:length(nr)
%     cov_m(1+I*(ii-1),1+I*(ii-1))=Pvari;
%     cov_m(2+I*(ii-1),2+I*(ii-1))=Svari;
%     cov_m(3+I*(ii-1),3+I*(ii-1))=Rvari;
%     
%     cov_m(2+I*(ii-1),1+I*(ii-1))=0; cov_m(1+I*(ii-1),2+I*(ii-1))=0;
%     cov_m(3+I*(ii-1),1+I*(ii-1))=0;  cov_m(1+I*(ii-1),3+I*(ii-1))=0; 
%     cov_m(3+I*(ii-1),2+I*(ii-1))=0;  cov_m(2+I*(ii-1),3+I*(ii-1))=0; 
% end


G_new_alpha=G_new(:,1:11);
G_alpha_vec = reshape(G_new_alpha,1,[])';

d_vec = reshape(d_new,1,[])';

obs_d_gradient_matrix = reshape(obs_d_gradient,140,11);

partt1 = cov_m*G_new_alpha';
partt2 = inv(G_new_alpha*cov_m*G_new_alpha'+cov_d);
partt3  = obs_d_gradient_matrix-G_new_alpha.*m_est_1D(:,1);
partt3=zeros(140,11);

mean_posterior = m_est_1D_alpha+partt1*partt2.*partt3';
mean_posterior=mean_posterior';

m_posterior = inv(G_new_alpha'*inv(cov_d)*G_new_alpha+inv(cov_m));

conf_area = diag(m_posterior);

%mean_posterior = m_est_new+cov_m*G_new'*inv(G_new*cov_m*G_new'+cov_d)*(d_new-G_new*m_est_new); % eq 3.37 tarantola
size(mean_posterior)

%% Endre fra mean_posterior til mean_posterior_vel

%mean_posterior_vel = zeros(length(time),3);
m_est_1D_alpha = m_est_1D(:,1);
    
m_inv_alpha_1D=zeros(length(time),1);

for i = 2:length(nr)
    m_inv_alpha_1D(1:nr(1)) = vpt(1);
    m_inv_alpha_1D(nr(i-1)+1:nr(i)) = m_est_1D_alpha(nr(i))*vpt(i-1)+vpt(i-1);
    m_inv_alpha_1D(nr(end):end) = m_est_1D_alpha(nr(end)+1)*vpt(end-1)+vpt(end-1);
end


mean_posterior_vel_test = zeros(length(time),1);

for i = 2:length(nr)
    mean_posterior_vel_test(1:nr(1)) = vpt(1);
    mean_posterior_vel_test(nr(i-1)+1:nr(i)) = mean_posterior(nr(i))*vpt(i-1)+vpt(i-1);
    mean_posterior_vel_test(nr(end):end) = mean_posterior(nr(end)+1)*vpt(end-1)+vpt(end-1);
end


%% Confidence area
I=3;
cov_d = eye(length(d))*max(max(Rd));
cov_m=eye(Num_parameters); %covariance model
for ii=1:length(nr)
    cov_m(1+I*(ii-1),1+I*(ii-1))=Pvari;
    cov_m(2+I*(ii-1),2+I*(ii-1))=Svari;
    cov_m(3+I*(ii-1),3+I*(ii-1))=Rvari;
    
    cov_m(2+I*(ii-1),1+I*(ii-1))=0; cov_m(1+I*(ii-1),2+I*(ii-1))=0;
    cov_m(3+I*(ii-1),1+I*(ii-1))=0;  cov_m(1+I*(ii-1),3+I*(ii-1))=0; 
    cov_m(3+I*(ii-1),2+I*(ii-1))=0;  cov_m(2+I*(ii-1),3+I*(ii-1))=0; 
end

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

%%

hold all
p1 = patch([vp1D_max fliplr(vp1D_min)], [time fliplr(time)],[0.2,0.2,0.2]); set(p1, 'facealpha',.5)
%p1 = patch([vp1D_max fliplr(vp1D_min)], [time fliplr(time)],'r'); set(p1, 'facealpha',.5)
stairs(mean_posterior_vel_test,time,'r','Linewidth',3.5)
stairs(vp1D,time,'k','Linewidth',3.5),title('P-wave velocity'),grid on
hold off
xlabel('Velocity (m/s)'),ylabel('TWT (s)')
set(gca,'Ydir','reverse'),set(gca,'FontSize',10),set(gca,'Linewidth',2)
%axis([x_min_vp,x_max_vp,y_min,y_max]),set(gca,'FontSize',35)

