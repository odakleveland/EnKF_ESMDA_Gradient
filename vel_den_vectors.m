function [vp1D,vs1D,rho1D] = vel_den_vectors(time,nr,vpt,vst,rhot)
%function that creates vectors for p, s and density plot

nnr=[0,nr,length(time)];
for gg=1:length(nnr)-1
    vp1D(1+nnr(gg):nnr(gg+1))=vpt(gg);
    vs1D(1+nnr(gg):nnr(gg+1))=vst(gg);
    rho1D(1+nnr(gg):nnr(gg+1))=rhot(gg);
end

end