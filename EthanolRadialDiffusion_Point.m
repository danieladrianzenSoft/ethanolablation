function EthanolRadialDiffusion_Point

clear
clc

[D_S,Rinit,Rtot,phi,C_0,t] = initializeParameters();
tspan = [min(t),max(t)];

numr = 800;
ns = numr;

r = linspace(0,Rtot,numr);
nRinit = find(r>=Rinit,1);
bound = find(r>=0.5,1); %Tumor boundary at 0.5cm radius (1cm diameter)

  
IC = zeros(1,numr);
IC(1:nRinit) = C_0;

S = calcsparsity(r,numr);
S = sparse(S);


opts1 = odeset('Stats','on','JPattern',S);

[t2,C] = ode15s(@(t2,C) dcdt(t,C,r,Rinit,Rtot,D_S,phi,ns,C_0), tspan, IC,opts1);

C = C';

figure()
tvec1=[find(t2>=0.5*60,1) find(t2>=10*60,1) find(t2>=30*60,1)];
for i=1:length(tvec1)
plot(r,C(:,tvec1(i)),'LineWidth',2)
hold all
end

%  ylimmin = 0;
%  ylimmax = 1;
%  axis([0 x(end) ylimmin ylimmax])
%  ylims = linspace(ylimmin, ylimmax, 500);
%  x1lim = h_S*(ylims./ylims);
%  hold on
%  plot(x1lim,ylims,'k--')
 ylabel('Concentration (normalized)','FontSize',18,'FontWeight','Bold')
 xlabel('Radius (cm)','FontSize',18,'FontWeight','Bold')
 %title('Normalized ethanol concentration vs. radius at different times','FontSize',18,'FontWeight','Bold')
 legend('0.5mins','10mins','30mins')
 XT = get(gca,'XTick');
 set(gca,'FontSize',16)
 YT = get(gca,'YTick');
 set(gca,'FontSize',16)
 xlim([0,0.8])
 
MassInside = trapz(r(1:bound),C(1:bound,:),1);
MassOutside = trapz(r(bound+1:numr),C(bound+1:numr,:),1);
%MassInside = trapz(C(1:bound,:),1);
%MassOutside = trapz(C(bound+1:end,:),1);
TotMass = MassInside+MassOutside;
TotMassPredic = trapz(r,C,1);

% figure()
% plot(t2/60,MassInside,'b','LineWidth',2)
% hold on
% plot(t2/60,MassOutside,'r','LineWidth',2)
% plot(t2/60,TotMass,'k','LineWidth',2)
% plot(t2/60,TotMassPredic,'g--','LineWidth',2)
% 
% legend('Inside Tumor','Outside Tumor','Total Drug Mass','Total Predicted')
 
%  figure()
%  plot(t2/60,C(bound,:),'LineWidth',2)
%  xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
%  ylabel('Ethanol concentration (normalized)','FontSize',18,'FontWeight','Bold')
%  title('Ethanol concentration at tumor boundary r = 0.5cm vs. time','FontSize',18,'FontWeight','Bold')
%  XT = get(gca,'XTick');
%  set(gca,'FontSize',16)
%  YT = get(gca,'YTick');
%  set(gca,'FontSize',16)

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

     X = linspace(min(data(:,1)),max(data(:,1)),150);
     Y = linspace(min(data(:,2)),max(data(:,2)),150);

     [XX,YY]=meshgrid(X,Y);
     
     figure()

     F=scatteredInterpolant(data(:,1),data(:,2),data(:,3));
     contourf(XX,YY,F(XX,YY),100,'LineColor','none')
     xlim([-1.2,1.2])
     ylim([-1.2,1.2])
     hold on
     colormap bone
     colorbar
     plot(tumorX,tumorY,'w--','LineWidth',2)
     xlabel('X (cm)','FontSize',18,'FontWeight','Bold')
     ylabel('Y (cm)','FontSize',18,'FontWeight','Bold')
     title(sprintf('t = %.0f mins',floor(t2(tvec1(j))/60)),'FontSize',18,'FontWeight','Bold')
     XT = get(gca,'XTick');
     set(gca,'FontSize',16)
     YT = get(gca,'YTick');
     set(gca,'FontSize',16)

     
 end

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

function [D_S,Rinit,Rtot,phi,C_0,t] = initializeParameters()

%%%%%% PARAMETERS %%%%%%

% D_OH = 1*10^(-5); %cm^2/s, from "Molecular Dynamics simulation of ethanol/water
% %mixtures for structure and diffusion properties, 2005.
D_S = 7.5*10^(-6); %(cm^2/s)
C_0 = 1;
phi = 1;

%%%GEOMETRY:

Rtot = 1.5; %cm - Limit of our analysis
Rinit = 0.25/2; %cm - Initial ethanol distribution
t = linspace(0,2*60*60,1001);

end



function DCDT = dcdt(t,C,r,Rinit,Rtot,D_S,phi,ns,C_0)

%NOTE: TO UNDERSTAND THESE EQUATIONS, LOOK UP FOR DIFFERENCE EQUATIONS: CENTRAL, BACKWARD
%AND FORWARD.
dr = r(2)-r(1);
dCdt = zeros(1,length(r));
dCdr = zeros(1,length(r));
d2Cdr2 = zeros(1,length(r));

for i = 1 %y direction BC zero flux (dcdx=0)
    dCdt(i) = 6*D_S.*(C(i+1)-C(i))./(dr.^2);
end
for i = 2:(ns-1) %in stroma
    %dCdr(i) = (C(i+1)-C(i))./dr;
    %dCdr(i) = (C(i)-C(i-1))./dr;
    dCdr(i) = (C(i+1)-C(i-1))./(2*dr);
    d2Cdr2(i) = (C(i+1)-2.*C(i)+C(i-1))./(dr.^2);
    dCdt(i) = ((2*D_S/r(i)).*dCdr(i))+(D_S*d2Cdr2(i));
end
for i = ns %right before interface - stroma
    dCdt(i) = 6*D_S*(C(i-1)-C(i))./(dr.^2);
    %dCdt(i) = 0;
end

DCDT = dCdt';

end

function S = calcsparsity(r,numr)
    totalSize = numr;
    S = ones(totalSize,totalSize);
end