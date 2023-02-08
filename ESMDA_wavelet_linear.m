function Xp=ESMDA_wavelet_linear(X,w,a,obs,R,Num_ensembles,Num_parameters,Na,bvpt,bvst,brhot,angl,time,nr)

%ES-MDA updates the ensembles X. From this the posterior covariance and mean
%is found.

%Input:
%X     - Prior state vector
%w     - Wavelet in time domain
%a     - Inflation coefficient
%obs   - observed data
%R     - measurment covariance matrix
%N     - number of ensemble members
%M     - number of parameters
%Na    - number of assimilations
%bvpt  - background p-wave vel
%bvst  - background s-wave vel
%brhot - background density
%angl  - angle
%time  - time
%nr    - number of interfaces

%Output:
%Xp    - posterior state vector


nt=length(time);                   %Length of the time/layer vector
I=3;                               %Number of parameters per layer
Num_obs=length(obs);                     %Total number of observations

Xp=X;
%Loop over all assimilations
%rng(13613)                                     %Seed number, not needed
for hh=1:Na    
    hh
    %Replicate observation vector
    
    drand=sqrt((a(hh))*R)*randn(Num_obs,Num_ensembles);              %Samples from a normal distribution ~ N(0,diag(R))
    d=zeros(Num_obs,Num_ensembles);
    for jj=1:Num_obs
        d(jj,:)=obs(jj)+drand(jj,:);               %Replicate observed data,adding noise twice
    end

    
    %Forward modeling operator. Loop over all interfaces, angles and ensembles
    %Forward modelling operator
    Hx=zeros(length(angl)*(nt),Num_ensembles); 
    pp=0;
    for jj=1:length(angl)
        Rpp1=zeros(length(time),Num_ensembles);
        for ii=1:Num_ensembles
            for ll=1:length(nr)
                %Reflection coefficient
                Rpp1(nr(ll),ii)=linear_zoeppritz_contrast(angl(jj),bvpt(ll),Xp(1+I*(ll-1),ii)-bvpt(ll),bvst(ll),...
                    Xp(2+I*(ll-1),ii)-bvst(ll),brhot(ll),Xp(3+I*(ll-1),ii)-brhot(ll));
                
            end
            %Convolution with the wavelet, check the wavelet
            uz_l1=conv(Rpp1(:,ii),w,'same');  %check conv function            
            uz=uz_l1(1:nt); %might not need this
            
            %Forward modeled data, change for ray born
            Hx((jj-1)*length(time)+1:jj*length(time),ii)=uz; 
        end
    end
    
    Hx=real(Hx);

    %See equations from the ES-MDA algorithm for more information.
    E=1/Num_ensembles*sum(Xp');                 %Expectation found from Xp
    A=zeros(Num_parameters,Num_ensembles);
    for kk=1:Num_parameters
        A(kk,:)=Xp(kk,:)-E(kk);     %Perturbation matrix A
    end
    
    Eh=1/Num_ensembles*sum(Hx');                %Expectation from data Hx
    HA=zeros(Num_obs,Num_ensembles);
    for ll=1:Num_obs
        HA(ll,:)=Hx(ll,:)-Eh(ll);   %Perturbation matrix X
    end

    P=1/(Num_ensembles-1)*HA*HA'+R*(a(hh));             %P matrix, from equation 4.30
    
    Xp=Xp+1/(Num_ensembles-1)*A*HA'*inv(P)*(d-Hx);      %Update ensembles
end

end


