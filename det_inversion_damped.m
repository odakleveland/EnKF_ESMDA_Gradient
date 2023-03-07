function [vp_inv, vs_inv, rho_inv,m_inv,m_est] = det_inversion_damped(G,d,nr,vpt,vst,rhot,num_alpha) 

I = eye(size(G'*G));

n_alpha = length(I);
alpha = logspace(-2,0.01,n_alpha);

m_est = inv(G'*G+alpha(num_alpha).^2*I)*(G'*d);

m_inv = zeros(length(nr),3);

%create size for m_inv
for i = 1:length(nr)
    m_inv(i,1:3) = [m_est(3*(i-1)+1)*vpt(i)+vpt(i),m_est(3*(i-1)+2)*vst(i)+vst(i),m_est(3*i)*rhot(i)+rhot(i)];
end

%create vectors for velocities after inversion
vp_inv = [2300, m_inv(:,1)'];
vs_inv = [1170, m_inv(:,2)'];
rho_inv = [2146, m_inv(:,3)'];

m_inv = [vp_inv',vs_inv',rho_inv']

end