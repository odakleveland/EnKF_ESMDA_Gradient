function [G,d] = lin_zoeppritz(vpt,vst,rhot,theta,nr)

delta_alpha = diff(vpt);delta_beta = diff(vst);delta_rho = diff(rhot);

G = zeros(length(theta)*length(nr),3*length(nr));

for i = 1:length(nr)
    for j = 1:length(theta)
        a_alpha(j,i) = 0.5*(1+(tand(theta(j)).^2));
        a_beta(j,i) = -4*(vst(i).^2)/(vpt(i).^2)*sind(theta(j)).^2;
        a_rho(j,i) = 0.5*(1-4*(vst(i).^2)/(vpt(i).^2)*sind(theta(j)).^2);
        d(j,i) = a_alpha(j)*(delta_alpha(i)/vpt(i))+a_beta(j)*(delta_beta(i)/vst(i))+a_rho(j)*(delta_rho(i)/rhot(i));
    end  
end

a_alpha = reshape(a_alpha,[],1);
a_beta = reshape(a_beta,[],1);
a_rho = reshape(a_rho,[],1);
d = reshape(d,[],1);

for i = 1:length(nr)
    G(((i-1)*11)+1:i*11,(3*(i-1))+1:3*i) = [a_alpha((i-1)*11+1:i*11),a_beta((i-1)*11+1:i*11),a_rho((i-1)*11+1:i*11)];
end

end