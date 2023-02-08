function Rpp = linear_zoeppritz_contrast(theta,vp1,vp2,vs1,vs2,rho1,rho2)

% Input:
% theta - angles
% vp1 - p-wave velocity, first layer
% vp2 - p-wave velocity, pertubation in second layer wrt first layer
% vs1 - s-wave velocity, first layer
% vs2 - s-wave velocity, pertubation in second layer wrt first layer
% rho1 - density, first layer
% rho2 - density, pertubation in second layer wrt first layer
% 
% Output:
% Rpp - reflection coefficients

vp22 = vp2 + vp1;
vs22 = vs2 + vs1;
rho22 = rho2 + rho1;

alpha1=(vp1+vp22)/2;
beta1=(vs1+vs22)/2;

%Linear zoeppritz:
a_alpha = 0.5*(1+(tand(theta).^2));
a_beta = -4*(beta1^2/alpha1^2)*(sind(theta).^2);
a_rho = 0.5*(1-4*(beta1^2/alpha1^2)*(sind(theta).^2));

delta_alpha = vp22-vp1;
delta_beta = vs22-vs1;
delta_rho = rho22-rho1;

Rpp = a_alpha*(delta_alpha/alpha1) + a_beta*(delta_beta/beta1) + a_rho*(delta_rho/rho1);
end