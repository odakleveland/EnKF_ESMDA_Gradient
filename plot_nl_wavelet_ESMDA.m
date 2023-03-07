function [pvp95min,pvp95max,pvs95min,pvs95max,pr95min,pr95max,x2,zz,vpEp,vsEp,rEp]= plot_nl_wavelet_ESMDA(X,Xp,vp,vs,rho,time,nr,I,distr)

%Input
%X    - prior state vector
%Xp   - posterior state vector
%vp   - true p-wave velocity
%vs   - true s-wave velocity
%rho  - true density
%time - time
%nr   - number of interfaces
%I    - parameters per layer
%distr- 1=gaussian dist.  2=uniform dist.

%No output


%Expectence and variance for prior and posterior
if distr==1
    E=mean(X'); %gaussian dist.
    Var=std(X');
end
if distr==2
    b=max(X'); a=min(X');
    E=(b-a)/2+a;
    Var=b-E;             %uniform dist.
end
% [Ep,Varp]=normfit(Xp');
Ep=mean(Xp');
Varp=std(Xp');

%Expectence and variance in the first layer is not known
Epvp(1)=vp(1);      Varpp(1)=0;     Evp(1)=E(1);       Varvp(1)=0; 
pmin95vp(1)=vp(1);  pmax95vp(1)=vp(1);  min95vp(1)=vp(1);   max95vp(1)=vp(1);

Epvs(1)=vs(1);      Varps(1)=0;     Evs(1)=E(2);       Varvs(1)=0;
pmin95vs(1)=vs(1);  pmax95vs(1)=vs(1);  min95vs(1)=vs(1);   max95vs(1)=vs(1);

Epr(1)=rho(1);      Varpr(1)=0;    Evr(1)=E(3);      Varvr(1)=0;
pmin95r(1)=rho(1);  pmax95r(1)=rho(1);  min95r(1)=rho(1);   max95r(1)=rho(1);

%Expectence, variance and 95% confidence interval in all layers below layer one
conf = 1.960; %conf = 1.960;, conf = 1;
for ii=1:length(vp)-1
    Epvp(ii+1)=Ep(1+I*(ii-1));    Varpp(ii+1)=Varp(1+I*(ii-1));     %P-vel
    Epvs(ii+1)=Ep(2+I*(ii-1));    Varps(ii+1)=Varp(2+I*(ii-1));     %S-vel
    Epr(ii+1)=Ep(3+I*(ii-1));     Varpr(ii+1)=Varp(3+I*(ii-1));     %Density
    
    Evp(ii+1)=E(1+I*(ii-1));    
    Evs(ii+1)=E(2+I*(ii-1));    
    Evr(ii+1)=E(3+I*(ii-1));    
    
    Varvp(ii+1)=Var(1+I*(ii-1));        %P-vel
    Varvs(ii+1)=Var(2+I*(ii-1));        %S-vel
    Varvr(ii+1)=Var(3+I*(ii-1));        %Density
    
    pmin95vp(ii+1)=Epvp(ii+1)-conf*Varpp(ii+1)/1;pmax95vp(ii+1)=Epvp(ii+1)+conf*Varpp(ii+1)/1;%P-vel
    pmin95vs(ii+1)=Epvs(ii+1)-conf*Varps(ii+1)/1;pmax95vs(ii+1)=Epvs(ii+1)+conf*Varps(ii+1)/1;%S-vel
    pmin95r(ii+1)=Epr(ii+1)-conf*Varpr(ii+1)/1;pmax95r(ii+1)=Epr(ii+1)+conf*Varpr(ii+1)/1;    %Density
    
    if distr==1
    min95vp(ii+1)=Evp(ii+1)-conf*Varvp(ii+1)/1;max95vp(ii+1)=Evp(ii+1)+conf*Varvp(ii+1)/1;    %P-vel
    min95vs(ii+1)=Evs(ii+1)-conf*Varvs(ii+1)/1;max95vs(ii+1)=Evs(ii+1)+conf*Varvs(ii+1)/1;    %S-vel
    min95r(ii+1)=Evr(ii+1)-conf*Varvr(ii+1)/1;max95r(ii+1)=Evr(ii+1)+conf*Varvr(ii+1)/1;      %Density
    end
    if distr==2
    min95vp(ii+1)=Evp(ii+1)-0.95*Varvp(ii+1)/1;max95vp(ii+1)=Evp(ii+1)+0.95*Varvp(ii+1)/1;    %P-vel
    min95vs(ii+1)=Evs(ii+1)-0.95*Varvs(ii+1)/1;max95vs(ii+1)=Evs(ii+1)+0.95*Varvs(ii+1)/1;    %S-vel
    min95r(ii+1)=Evr(ii+1)-0.95*Varvr(ii+1)/1;max95r(ii+1)=Evr(ii+1)+0.95*Varvr(ii+1)/1;      %Density
    end
end

%Mean std in all layers
[mean(Varpp(2:end)),mean(Varps(2:end)),mean(Varpr(2:end))]

%rms mean in all layers
[rms(Epvp(2:end)-vp(2:end)),rms(Epvs(2:end)-vs(2:end)),rms(Epr(2:end)-rho(2:end))]


%smoothing of plot
min95vp(2:end)=min95vp(2:2);
max95vp(2:end)=max95vp(2:2);
Evp(2:end)=Evp(2);
min95vs(2:end)=min95vs(2:2);
max95vs(2:end)=max95vs(2:2);
Evs(2:end)=Evs(2);
min95r(2:end)=min95r(2:2);
max95r(2:end)=max95r(2:2);
Evr(2:end)=Evr(2);
% Epvp=smooth(Epvp,7);
nnr=[0,nr-1,length(time)-1];

%For fill plot
for hh=1:length(min95vp)
    vp95min(1+nnr(hh):nnr(hh+1))=min95vp(hh);
    pvp95min(1+nnr(hh):nnr(hh+1))=pmin95vp(hh);
    vp95max(1+nnr(hh):nnr(hh+1))=max95vp(hh);
    pvp95max(1+nnr(hh):nnr(hh+1))=pmax95vp(hh);
    vpEp(1+nnr(hh):nnr(hh+1))=Epvp(hh);
    vpE(1+nnr(hh):nnr(hh+1))=Evp(hh);
    vpp(1+nnr(hh):nnr(hh+1))=vp(hh);
    
    vs95min(1+nnr(hh):nnr(hh+1))=min95vs(hh);
    pvs95min(1+nnr(hh):nnr(hh+1))=pmin95vs(hh);
    vs95max(1+nnr(hh):nnr(hh+1))=max95vs(hh);
    pvs95max(1+nnr(hh):nnr(hh+1))=pmax95vs(hh);
    vsEp(1+nnr(hh):nnr(hh+1))=Epvs(hh);
    vsE(1+nnr(hh):nnr(hh+1))=Evs(hh);
    vsp(1+nnr(hh):nnr(hh+1))=vs(hh);
    
    r95min(1+nnr(hh):nnr(hh+1))=min95r(hh);
    pr95min(1+nnr(hh):nnr(hh+1))=pmin95r(hh);
    r95max(1+nnr(hh):nnr(hh+1))=max95r(hh);
    pr95max(1+nnr(hh):nnr(hh+1))=pmax95r(hh);
    rEp(1+nnr(hh):nnr(hh+1))=Epr(hh);
    rE(1+nnr(hh):nnr(hh+1))=Evr(hh);
    rp(1+nnr(hh):nnr(hh+1))=rho(hh);
end
zz=linspace(time(1),time(end),length(vp95min));


%Figures

x_min=min(min95vp)-200*1.96;
x_max=max(max95vp)+300*1.96;
y_min=time(nr(1)+1);
y_max=max(time);

%figure(112)
subplot(1,3,1)
axis([x_min,x_max,y_min,y_max]),hold on,set(gca,'FontSize',35)
set(gca,'Linewidth',4)
title('P-wave velocity')
xlabel('Velocity (m/s)'),ylabel('TWT (s)')
set(gca,'Ydir','reverse')
x2=[zz,fliplr(zz)];
inBetween1=[vp95min,fliplr(vp95max)];
inBetween2=[pvp95min,fliplr(pvp95max)];
h1=fill(inBetween1,x2,[0.8,0.8,0.8]); set(h1,'facealpha',.5)
h2=fill(inBetween2,x2,[0.2,0.2,0.2]); set(h2,'facealpha',.55)
stairs(vpE,zz,'b','Linewidth',3.5)
stairs(vpEp,zz,'k','Linewidth',3.5)
stairs(vpp,zz,'r','Linewidth',3.5),grid on
%legend('95% pr. conf.','95% post. conf.','Prior mean','Post. mean','True model')
hold off

x_min=min(min95vs)-200*1.96;
x_max=max(max95vs)+200*1.96;
y_min=time(nr(1)+1);
y_max=max(time);

%figure(113)
subplot(1,3,2)
axis([x_min,x_max,y_min,y_max]),hold on,set(gca,'FontSize',35)
set(gca,'Linewidth',4)
title('S-wave velocity')
xlabel('Velocity (m/s)'),ylabel('TWT (s)')
set(gca,'Ydir','reverse')
x2=[zz,fliplr(zz)];
inBetween1=[vs95min,fliplr(vs95max)];
inBetween2=[pvs95min,fliplr(pvs95max)];
h3=fill(inBetween1,x2,[0.8,0.8,0.8]); set(h3,'facealpha',.5)
h4=fill(inBetween2,x2,[0.2,0.2,0.2]); set(h4,'facealpha',.55)
stairs(vsE,zz,'b','Linewidth',3.5)
stairs(vsEp,zz,'k','Linewidth',3.5)
stairs(vsp,zz,'r','Linewidth',3.5),grid on
%legend('95% pre. prediction','95% post. prediction','Prior','Posterior','True model')
hold off


x_min=min(min95r)-200*1.96;
x_max=max(max95r)+200*1.96;
y_min=time(nr(1)+1);
y_max=max(time);

%figure(114)
subplot(1,3,3)
axis([x_min,x_max,y_min,y_max]),hold on,set(gca,'FontSize',35)
set(gca,'Linewidth',4)
title('Density')
xlabel('Density (kg/m^3)'),ylabel('TWT (s)')
set(gca,'Ydir','reverse')
x2=[zz,fliplr(zz)];
inBetween1=[r95min,fliplr(r95max)];
inBetween2=[pr95min,fliplr(pr95max)];
h5=fill(inBetween1,x2,[0.8,0.8,0.8]); set(h5,'facealpha',.5)
h6=fill(inBetween2,x2,[0.2,0.2,0.2]); set(h6,'facealpha',.55) %post conf
stairs(rE,zz,'b','Linewidth',3.5)
stairs(rEp,zz,'k','Linewidth',3.5) %post mean
stairs(rp,zz,'r','Linewidth',3.5),grid on %true model
hb=legend('95% prior confidence area','95% posterior confidence area','Prior mean','Posterior mean','True model')
hold off

save hb

end