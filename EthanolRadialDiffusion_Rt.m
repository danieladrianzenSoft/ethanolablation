function EthanolRadialDiffusion_Rt

clear
clc

[D_S,D_C,Rfin,Rtot,R0,phi,C_0,ti,tau,totVol] = initializeParameters();

%t = linspace(0,2*60*60,801);
%tspan = [min(ti),max(ti)];

numr = 1601;
ns = numr;
r = linspace(0,Rtot,numr);

t2 = linspace(0,2*60*60,1201);
Retcell = (Rfin-R0)*(1-exp(-t2/tau))+R0;
nRetcell = zeros(1,length(t2));

%t2 = zeros(1,length(Ttemp));
%Retcell = interp1(Ttemp,Retcell_temp,t);

bound = find(r>=0.5,1); %Tumor boundary at 0.5cm radius (1cm diameter)

IC = zeros(1,numr); %Initial 'initial condition'
%IC(1:nRetcell) = C_0;

S = calcsparsity(r,numr);
S = sparse(S);
 
%opts1 = odeset('Stats','on','JPattern',S);
opts1 = odeset('JPattern',S);

C = zeros(length(t2),length(r));
te = 0;

tempRetcell = (Rfin-Retcell)/Rfin;
index = find(tempRetcell<=0.0005,1); %%Index in time when radius stops growing.

for jj = 1:index

    nRetcell(jj) = find(r>=Retcell(jj),1);
    if jj == 1
        t = [0,ti];
        %nRetcell(jj) = 2;
    else
        t = te+ti;
        %nRetcell(jj) = find(r>=Retcell(jj),1);
        IC(1:nRetcell(jj)) = C_0;
    end

    [t,Ct] = ode15s(@(t,Ct) dcdt(t,Ct,r,Rfin,Rtot,D_S,D_C,phi,ns,tau,nRetcell(jj),C_0), t, IC,opts1);
    te = t(end);
    IC = Ct(end,:);
    
    %t2 = [t2;t];
    C(jj,:) = mean(Ct,1);

end

remainTime = t2(index+1:end);
nRetcell(index+1:end) = nRetcell(index);
[t,Ct] = ode15s(@(t,Ct) dcdt(remainTime,Ct,r,Rfin,Rtot,D_S,D_C,phi,ns,tau,nRetcell(index),C_0), remainTime, IC,opts1);
C(index+1:end,:) = Ct;
C = C';

% figure()
% plot(t2(1:index)/60,r(nRetcell(1:index)),'LineWidth',2)
% hold on
% plot(t2(1:index)/60,Retcell(1:index),'g--','LineWidth',2)
% ylabel('Radius of Ethylcellulose - Ethanol cavity','FontSize',18,'FontWeight','Bold')
% xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
% title('Radius of cavity over time during injection','FontSize',18,'FontWeight','Bold')
% XT = get(gca,'XTick');
% set(gca,'FontSize',16)
% YT = get(gca,'YTick');
% set(gca,'FontSize',16)
% legend('Model','Defined function')
%  
% figure()
% plot(t2(1:index)/60,Retcell(1:index),'LineWidth',2)
% ylabel('Radius of Ethylcellulose - Ethanol cavity','FontSize',18,'FontWeight','Bold')
% xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
% title('Radius of cavity over time during injection','FontSize',18,'FontWeight','Bold')
% XT = get(gca,'XTick');
% set(gca,'FontSize',16)
% YT = get(gca,'YTick');
% set(gca,'FontSize',16)
%  

MassInside = 4*pi*trapz(r(1:bound),(r(1:bound)'.^2).*C(1:bound,:),1);
MassOutside = 4*pi*trapz(r(bound+1:numr),(r(bound+1:numr)'.^2).*C(bound+1:numr,:),1);
TotMass = MassInside+MassOutside;
TotMassPredic = (C_0*totVol)*ones(1,length(t2));

figure()
plot(t2/60,MassInside,'LineWidth',2)
hold on
plot(t2/60,MassOutside,'r','LineWidth',2)
plot(t2/60,TotMass,'k--','LineWidth',2)
ylabel('Ethanol Mass (mg)','FontSize',18,'FontWeight','Bold')
xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
title('Total Mass of Ethanol Inside and Outside Tumor','FontSize',18,'FontWeight','Bold')
XT = get(gca,'XTick');
set(gca,'FontSize',16)
YT = get(gca,'YTick');
set(gca,'FontSize',16)
legend('Inside tumor','Outside tumor','Total')
% plot(t2/60,TotMassPredic,'g--','LineWidth',2)
% ylim([0,0.1])

figure()

tvec1=[find(t2>=0.5*60,1) find(t2>=20*60,1) find(t2>=30*60,1)];

cm = colormap('copper');
colordivision = floor(size(cm,1)/length(tvec1));

for i=1:length(tvec1)

plot(r,C(:,tvec1(i)),'color',cm(i*colordivision,:),'LineWidth',2)
hold all

end

 ylabel('Concentration (normalized)','FontSize',18,'FontWeight','Bold')
 xlabel('Radius (cm)','FontSize',18,'FontWeight','Bold')
 %title('Normalized ethanol concentration vs. radius at different times','FontSize',18,'FontWeight','Bold')
 legend('0.5mins','20mins','30mins')
 XT = get(gca,'XTick');
 set(gca,'FontSize',16)
 YT = get(gca,'YTick');
 set(gca,'FontSize',16)
 xlim([0,0.8])


%  figure()
% 
%  plot(t2/60,C(bound,:),'LineWidth',2)
%  xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
%  ylabel('Ethanol concentration (normalized)','FontSize',18,'FontWeight','Bold')
%  title('Ethanol concentration at tumor boundary r = 0.5cm vs. time','FontSize',18,'FontWeight','Bold')
%  XT = get(gca,'XTick');
%  set(gca,'FontSize',16)
%  YT = get(gca,'YTick');
%  set(gca,'FontSize',16)
 
 makeSpherePlots(t2,tvec1,r,C,bound)


%  CTumor_avg = trapz(C(1:bound,:))/(ns-1);
%  figure()
%  plot(t2/60,CTumor_avg,'LineWidth',2)
%  xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
%  ylabel('Avg. Ethanol concentration (normalized)','FontSize',18,'FontWeight','Bold')
%  title('Average Ethanol concentration in tumor vs. time','FontSize',18,'FontWeight','Bold')
%  XT = get(gca,'XTick');
%  set(gca,'FontSize',16)
%  YT = get(gca,'YTick');
%  set(gca,'FontSize',16)
%  xlim([0,5])

% figure()
% plot(t2/60,PsurfaceC,'LineWidth',2)
% hold on
% plot(measurements(:,1),MsurfaceC,'x')
% %plot(TTemp,ConcTemp,'ro')
% legend('Model Predictions','Raman Measurements')
% xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
% ylabel('Ethanol concentration (normalized)','FontSize',18,'FontWeight','Bold')
% title('Ethanol concentration vs time on stromal surface','FontSize',18,'FontWeight','Bold')
% XT = get(gca,'XTick');
% set(gca,'FontSize',16)
% YT = get(gca,'YTick');
% set(gca,'FontSize',16)

end

 

function [D_S,D_C,Rfin,Rtot,R0,phi,C_0,ti,tau,totVol] = initializeParameters() 

%%%%%% PARAMETERS %%%%%%

% D_OH = 1*10^(-5); %cm^2/s, from "Molecular Dynamics simulation of ethanol/water
% %mixtures for structure and diffusion properties, 2005.

D_S = 7.5*10^(-6); %(cm^2/s)
D_C = D_S;
C_0 = 1;
phi = 1;
%%%GEOMETRY:

%t = linspace(0,2*60*60,1001);
ti = linspace(1,6,6);
Rtot = 2; %cm - Limit of our analysis
Rfin = 0.15; %cm - Ethanol distribution final radius after injection
R0 = 0.01; %Ethanol distribution initial radius after injection
flowRate = 1/(3600); %1ml/hr converted to ml/sec
totVol = 0.10; %mls (=100ul)
tau = (totVol/flowRate)/4; %time constant for increase of Retcell, divided 
%by 4 so that 98.2% of Rinit is achieved by the time injection ends.
%Remaining 1.8% expands in time as pressures equilibrate.

end

function DCDT = dcdt(t,C,r,Rinit,Rtot,D_S,D_C,phi,ns,tau,nRetcell,C_0)

%NOTE: TO UNDERSTAND THESE EQUATIONS, LOOK UP FOR DIFFERENCE EQUATIONS: CENTRAL, BACKWARD
%AND FORWARD.

dr = r(2)-r(1);
dCdt = zeros(1,length(r));
dCdr = zeros(1,length(r));
d2Cdr2 = zeros(1,length(r));

%C at cavity/stroma interface
C_intf_CSa = (D_C.*C(nRetcell)+(D_S.*C(nRetcell+1)))./((D_S.*phi)+D_C); %Right before interface
C_intf_CSb = (D_C.*C(nRetcell)+(D_S.*C(nRetcell+1)))./((D_C./phi)+D_S); %Right after interface

for i = 1 %y direction BC zero flux (dcdx=0)
    dCdt(i) = 6*D_C.*(C(i+1)-C(i))./(dr.^2);
end
for i = 2:nRetcell-1 %inside cavity
    dCdr(i) = (C(i+1)-C(i))./dr;
    d2Cdr2(i) = (C(i+1)-2.*C(i)+C(i-1))./(dr.^2);
    dCdt(i) = ((2*D_C/r(i)).*dCdr(i))+(D_C*d2Cdr2(i));
end
for i = nRetcell %right before interface - cavity/stroma
   dCdr(i) = (C_intf_CSa - C(i))./dr;
   d2Cdr2(i) = (C(i-1)-2.*C(i)+C_intf_CSa)./(dr.^2);
   dCdt(i) = ((2*D_C/r(i)).*dCdr(i))+(D_C*d2Cdr2(i));
end
for i = (nRetcell+1) %right after interface - cavity/stroma
    dCdr(i) = (C(i+1)-C(i))./dr;
    d2Cdr2(i) = ((C_intf_CSb)-2.*C(i)+C(i+1))./(dr.^2);
    dCdt(i) = ((2*D_S/r(i)).*dCdr(i))+(D_S*d2Cdr2(i));
end
for i = nRetcell+2:(ns-1) %in stroma
    dCdr(i) = (C(i+1)-C(i))./dr;
    d2Cdr2(i) = (C(i+1)-2.*C(i)+C(i-1))./(dr.^2);
    dCdt(i) = ((2*D_S/r(i)).*dCdr(i))+(D_S*d2Cdr2(i));
end
for i = ns %end of stroma, B.C.
    dCdt(i) = 6*D_S*(C(i-1)-C(i))./(dr.^2);
end

DCDT = dCdt';

end

function makeSpherePlots(t2,tvec1,r,C,bound)

 %CONVERTING TO CARTESIAN FOR PLOTTING
 
 theta = linspace(0,2*pi,200);
 x = zeros(1+(length(r)-1)*length(theta),1);
 y = zeros(1+(length(r)-1)*length(theta),1);
 z = zeros(1+(length(r)-1)*length(theta),1);
 
 tumorX = r(bound)*cos(theta);
 tumorY = r(bound)*sin(theta);
 
 for j = 1:length(tvec1)
     
     for i = 1:length(theta)
         x((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = r(2:end)*cos(theta(i));
         y((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = r(2:end)*sin(theta(i));
         z((length(r)-1)*(i-1)+2:(length(r)-1)*i+1) = C(2:end,tvec1(j));
     end

     z(1) = C(1,tvec1(j));

     data = [x,y,z];

     X = linspace(min(data(:,1)),max(data(:,1)),300);
     Y = linspace(min(data(:,2)),max(data(:,2)),300);

     [XX,YY]=meshgrid(X,Y);
     
     figure()

     F=scatteredInterpolant(data(:,1),data(:,2),data(:,3));
     contourf(XX,YY,F(XX,YY),250,'LineColor','none')
     xlim([-1.2,1.2])
     ylim([-1.2,1.2])
     hold on
     colormap parula
     colorbar
     caxis([0, 1]);
     plot(tumorX,tumorY,'w--','LineWidth',2)
     xlabel('X (cm)','FontSize',18,'FontWeight','Bold')
     ylabel('Y (cm)','FontSize',18,'FontWeight','Bold')
     title(sprintf('t = %.0f mins',floor(t2(tvec1(j))/60)),'FontSize',18,'FontWeight','Bold')
     XT = get(gca,'XTick');
     set(gca,'FontSize',16)
     YT = get(gca,'YTick');
     set(gca,'FontSize',16)
     
 end
end
 

function S = calcsparsity(r,numr)

    totalSize = numr;
    S = ones(totalSize,totalSize);

end

 
