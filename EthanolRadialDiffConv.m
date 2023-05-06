function EthanolRadialDiffConv

clear
clc

[D_S,Rinit,Rtot,phi,C_0,Vol,Q0,t2] = initializeParameters();
tspan = [min(t2),max(t2)];

numr = 1601;
ns = numr;

r = linspace(0,Rtot,numr);

nRinit = find(r>=Rinit,1);
bound = find(r>=0.5,1); %Tumor boundary at 0.5cm radius (1cm diameter)

tin_stop = find(t2>=Vol/Q0,1);
tpre = t2(1:tin_stop);
C = zeros(length(t2),length(r));

IC = zeros(1,numr);
IC(1:nRinit) = C_0;

S = calcsparsity(r,numr);
S = sparse(S);


opts1 = odeset('Stats','on','JPattern',S);

v = [0,Q0./(4*pi*r(2:end).^2)];
dvdr = [0,-Q0./(2*pi*r(2:end).^3)];
%v = zeros(1,length(r));
%dvdr = zeros(1,length(r));
[ttemp,Cpre] = ode15s(@(ttemp,Cpre) dcdt(tpre,Cpre,r,Rinit,Rtot,D_S,phi,ns,v,dvdr,C_0), tpre, IC,opts1);

avgVr = trapz(v)/(numr-1);
%avgVr = mean(v)
Pe = (avgVr*Rtot)/D_S;

C(1:tin_stop,:) = Cpre;

IC2 = C(tin_stop,:);
remaintime = t2(tin_stop+1:end);
v = zeros(1,length(r));
dvdr = zeros(1,length(r));
[ttemp2,Cpost] = ode15s(@(ttemp2,Cpost) dcdt(remaintime,Cpost,r,Rinit,Rtot,D_S,phi,ns,v,dvdr,C_0), remaintime, IC2,opts1);
C(tin_stop+1:end,:) = Cpost;

C = C';

figure()
cm = colormap('copper');

tvec1=[find(t2>=0*60,1)  find(t2>=1*60,1) find(t2>=2*60,1) find(t2>=4*60,1) find(t2>=6*60,1) find(t2>=10*60,1) find(t2>=15*60,1) find(t2>=30*60,1) find(t2>=60*60,1) find(t2>=2*60*60,1)];

%plot(r,C(:,tvec1(i)),'color',cm(i*6,:),'LineWidth',2)
%hold all
for i=1:length(tvec1)
plot(r,C(:,tvec1(i)),'color',cm(i*6,:),'LineWidth',2)
hold all
end

 ylabel('Ethanol Concentration (normalized)','FontSize',18,'FontWeight','Bold')
 xlabel('Radius (cm)','FontSize',18,'FontWeight','Bold')
 title('Normalized ethanol concentration vs. radius at different times','FontSize',18,'FontWeight','Bold')
 legend('0mins','1mins','2mins','4mins','6mins','10mins','15mins','30mins','60mins','120mins')
 XT = get(gca,'XTick');
 set(gca,'FontSize',16)
 YT = get(gca,'YTick');
 set(gca,'FontSize',16)
 xlim([0,0.8])
 ylim([0,1])
 
MassInside = 4*pi*trapz(r(1:bound),(r(1:bound)'.^2).*C(1:bound,:),1);
MassOutside = 4*pi*trapz(r(bound+1:numr),(r(bound+1:numr)'.^2).*C(bound+1:numr,:),1);
TotMass = MassInside+MassOutside;
TotMassPredic = (C_0*Vol)*ones(1,length(t2));

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
%plot(t2/60,TotMassPredic,'g--','LineWidth',2)
ylim([0,0.1])

legend('Inside Tumor','Outside Tumor','Total Injected')

% figure()
% plot(r,(Vr.*r)/D_S,'LineWidth',2)
%  xlabel('Radius (cm)','FontSize',18,'FontWeight','Bold')
%  ylabel('Peclet Number','FontSize',18,'FontWeight','Bold')
%  title('Peclet number vs radius','FontSize',18,'FontWeight','Bold')

 figure()
 plot(t2/60,C(bound,:),'LineWidth',2)
 xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
 ylabel('Ethanol concentration (normalized)','FontSize',18,'FontWeight','Bold')
 title('Ethanol concentration at tumor boundary r = 0.5cm vs. time','FontSize',18,'FontWeight','Bold')
 XT = get(gca,'XTick');
 set(gca,'FontSize',16)
 YT = get(gca,'YTick');
 set(gca,'FontSize',16)

[CM]=cbrewer('seq', 'Purples', 50, 'cubic');
makeSpherePlots(t2,tvec1,r,C,bound,CM)

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

end

function [D_S,Rinit,Rtot,phi,C_0,Vol,Q0,t2] = initializeParameters()

%%%%%% PARAMETERS %%%%%%

% D_OH = 1*10^(-5); %cm^2/s, from "Molecular Dynamics simulation of ethanol/water
% %mixtures for structure and diffusion properties, 2005.
D_S = 7.5*10^(-6); %(cm^2/s)
C_0 = 1;
phi = 1;
Vol = 0.1; %ml
Q0 = 1/3600; %1ml/hr to ml/s.

%%%GEOMETRY:

Rtot = 2; %cm - Limit of our analysis
%Rinit = 0.25/2; %cm - Initial ethanol distribution
Rinit = 0.01; %initial radius (27 gauge needle)
t2 = linspace(0,2*60*60,2401);

end

function DCDT = dcdt(t,C,r,Rinit,Rtot,D_S,phi,ns,v,dvdr,C_0)

%NOTE: TO UNDERSTAND THESE EQUATIONS, LOOK UP FOR DIFFERENCE EQUATIONS: CENTRAL, BACKWARD
%AND FORWARD.

dr = r(2)-r(1);
dCdt = zeros(1,length(r));
dCdr = zeros(1,length(r));
d2Cdr2 = zeros(1,length(r));

for i = 1 %y direction BC zero flux (dcdx=0, dvdr = 0)
   dCdt(i) = 6*D_S.*(C(i+1)-C(i))./(dr.^2);
end
for i = 2:(ns-1) %in tissue
    dCdr(i) = (C(i)-C(i-1))./dr;
    %dCdr(i) = (C(i+1)-C(i-1))./(2*dr);
    d2Cdr2(i) = (C(i+1)-2.*C(i)+C(i-1))./(dr.^2);
    dCdt(i) = ((2*D_S/r(i)).*dCdr(i))+(D_S*d2Cdr2(i))-((2/r(i)).*(v(i).*C(i)))-(C(i).*dvdr(i))-(v(i).*dCdr(i));
end
for i = ns %right before interface - stroma
    dCdt(i) = 6*D_S*(C(i-1)-C(i))./(dr.^2);
end

DCDT = dCdt';

end

function makeSpherePlots(t2,tvec1,r,C,bound,CM)

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
     xlim([-0.6,0.6])
     ylim([-0.6,0.6])
     hold on
     %colormap(flipud(CM));
     colormap('bone')
     colorbar;
     plot(tumorX,tumorY,'w--','LineWidth',2)
     xlabel('X (cm)','FontSize',24,'FontWeight','Bold')
     ylabel('Y (cm)','FontSize',24,'FontWeight','Bold')
     title(sprintf('Normalized Ethanol Concentration at %.0f mins',floor(t2(tvec1(j))/60)),'FontSize',18,'FontWeight','Bold')
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
 