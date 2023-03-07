function Rpp = reflection_coefficients(vp1D,vs1D,rho1D,theta)

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

end