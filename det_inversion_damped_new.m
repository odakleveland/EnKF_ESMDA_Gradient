function [vp_inv, vs_inv, rho_inv,m_inv,m_est] = det_inversion_damped_new(G,d,nr,vpt,vst,rhot,num_alpha) 

I = eye(size(G'*G));

n_alpha = length(I);
alpha = logspace(-2,0.01,n_alpha);

m_est = inv(G'*G+alpha(num_alpha).^2*I)*(G'*d);

m_inv = zeros(length(vpt),3);

for i = 1:length(vpt)
    m_inv(i,1:3) = [m_est(3*(i-1)+1)*vpt(i)+vpt(i),m_est(3*(i-1)+2)*vst(i)+vst(i),m_est(3*i)*rhot(i)+rhot(i)];
end

vp_inv = [m_inv(:,1)'];
vs_inv = [m_inv(:,2)'];
rho_inv = [m_inv(:,3)'];

m_inv = [vp_inv,vs_inv,rho_inv];

end