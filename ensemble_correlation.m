function [X]=ensemble_correlation(mu,Q,M,N,dist)

%Ensemble from a known pdf, assume Gaussian

%Input:
%E - expectation from prior
%Q - covariance from prior
%M - number of parameters
%N - number of ensembles
%dist - type of distribution, 1 is gaussian, 2 is uniform. 

%Output:
%X   - prior state vector

rng(1364)
if dist==1
    randnn=randn(M,N);              %Gaussian dist.
end
if dist==2
    randnn=rand(M,N)*2-1;           %Uniform
end
X=zeros(M,N);
R=chol(Q.^2);
X=R*randnn;
for kk=1:M
    X(kk,:)=X(kk,:)+mu(kk);   %Ensembles
end
end
