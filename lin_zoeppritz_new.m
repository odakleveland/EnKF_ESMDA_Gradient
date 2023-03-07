function [G,d,a_alpha,a_beta,a_rho] = lin_zoeppritz_new(vpt,vst,rhot,theta,nr,time)

delta_alpha = diff(vpt);delta_beta = diff(vst);delta_rho = diff(rhot);
delta_alpha(end+1)=0;delta_beta(end+1)=0;delta_rho(end+1)=0;

G = zeros(length(time),3*length(nr));

for i = 1:length(time)
    for j = 1:length(theta)
        a_alpha(i,j) = 0.5*(1+(tand(theta(j)).^2));
        a_beta(i,j) = -4*(vst(i).^2)/(vpt(i).^2)*sind(theta(j)).^2;
        a_rho(i,j) = 0.5*(1-4*(vst(i).^2)/(vpt(i).^2)*sind(theta(j)).^2);
        d(i,j) = a_alpha(1,j)*(delta_alpha(i)/vpt(i))+a_beta(1,j)*(delta_beta(i)/vst(i))+a_rho(1,j)*(delta_rho(i)/rhot(i));
    end  
end

G = [a_alpha,a_beta,a_rho];

end