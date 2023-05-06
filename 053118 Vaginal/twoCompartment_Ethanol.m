function twoCompartment_Ethanol

clear
clc

measurements = load('DataVsTime.csv');

[D_OH,D_S,k_E,h_OH,h_S,phi,C_0,t] = initializeParameters();

numx = 1000;
nOH = ceil((h_OH/(h_OH+h_S))*numx);
ns = numx-nOH;

x = [linspace(0,h_S,ns), linspace(h_S+(h_OH/nOH),h_S+h_OH,nOH)];
  
IC = zeros(1,nOH+ns);
IC(ns+1:numx) = C_0;

S = calcsparsity(x,numx,nOH,ns);
S = sparse(S);


opts1 = odeset('Stats','on','JPattern',S);

[t2,C] = ode15s(@(t2,C) dcdt(t,C,x,h_OH,h_S,D_OH,D_S,k_E,phi,nOH,ns,C_0), t, IC,opts1);

C = C';

%PsurfaceC = C(1,:)/mean(C(1,:));
PsurfaceC = C(1,:)/mean(C(1,:));

%measurements(:,1)
%round(t2/60)
indexes = zeros(length(measurements(:,1)),1);
for i = 1:length(measurements(:,1))
indexes(i) = find(round(t2/60)==measurements(i,1),1);
end
TTemp = round(t2(indexes,1)/60);
ConcTemp = PsurfaceC(indexes)'; %%fit

MsurfaceC = measurements(:,2)/mean(measurements(:,2));
MsurfaceCbar = mean(MsurfaceC);

SStot = sum((MsurfaceC-MsurfaceCbar).^2);
SSreg = sum((ConcTemp-MsurfaceCbar).^2);
SSres = sum((MsurfaceC-ConcTemp).^2);

R2 = 1-SSres/SStot

figure()
tvec1=[find(t>=5*60,1) find(t>=15*60,1) find(t>=30*60,1) find(t>=60*60,1) find(t>=2*60*60,1)];
for i=1:length(tvec1)
plot(x,C(:,tvec1(i)),'LineWidth',2)
hold all
end
 ylimmin = 0;
 ylimmax = 1;
 axis([0 x(end) ylimmin ylimmax])
 ylims = linspace(ylimmin, ylimmax, 500);
 x1lim = h_S*(ylims./ylims);
 hold on
 plot(x1lim,ylims,'k--')
 ylabel('Ethanol Concentration (normalized)','FontSize',18,'FontWeight','Bold')
 xlabel('Depth (cm)','FontSize',18,'FontWeight','Bold')
 title('Normalized ethanol vs. depth at different times','FontSize',18,'FontWeight','Bold')
 legend('5mins','15mins','30mins','60mins','120mins')
 XT = get(gca,'XTick');
 set(gca,'FontSize',16)
 YT = get(gca,'YTick');
 set(gca,'FontSize',16)

figure()
plot(t2/60,PsurfaceC,'LineWidth',2)
hold on
plot(measurements(:,1),MsurfaceC,'x')
%plot(TTemp,ConcTemp,'ro')
legend('Model Predictions','Raman Measurements')
xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
ylabel('Ethanol concentration (normalized)','FontSize',18,'FontWeight','Bold')
title('Ethanol concentration vs time on stromal surface','FontSize',18,'FontWeight','Bold')
XT = get(gca,'XTick');
set(gca,'FontSize',16)
YT = get(gca,'YTick');
set(gca,'FontSize',16)

end

function [D_OH,D_S,k_E,h_OH,h_S,phi,C_0,t] = initializeParameters()

%%%%%% PARAMETERS %%%%%%

D_OH = 1*10^(-5); %cm^2/s, from "Molecular Dynamics simulation of ethanol/water
%mixtures for structure and diffusion properties, 2005.
D_S = 5.3*10^(-6); %(cm^2/s) %FIT THIS
%D_S = 4.4*10^(-6);
C_0 = 1;
phi = 1;
k_E = 4.1/3600; %FIT THIS.
%k_E = 4/3600;
%%%GEOMETRY:

thickness = [0.8,0.75,0.89,0.8,0.73,0.75];
h_S = mean(thickness/10) %convert mm to cm %thickness tissue
h_OH = 1; %cm %thickness ethanol
t = linspace(0,81*60,1001);

end



function DCDT = dcdt(t,C,x,h_OH,h_S,D_OH,D_S,k_E,phi,nOH,ns,C_0)

%NOTE: TO UNDERSTAND THESE EQUATIONS, LOOK UP FOR DIFFERENCE EQUATIONS: CENTRAL, BACKWARD
%AND FORWARD.
dx = x(2)-x(1);
dCdt = zeros(1,length(x));
dCdx = zeros(1,length(x));
d2Cdx2 = zeros(1,length(x));

C_intfa = (D_S.*C(ns)+(D_OH.*C(ns+1)))./((D_OH.*phi)+D_S); %C at gel/tissue interface
C_intfb = (D_S.*C(ns)+(D_OH.*C(ns+1)))./((D_S./phi)+D_OH);

for i = 1 %y direction BC zero flux (dcdx=0)
    dCdt(i) = 2*D_S.*(C(i+1)-C(i))./(dx.^2)- k_E.*C(i);
end
for i = 2:(ns-1) %in stroma
    dCdx(i) = (C(i+1)-C(i))./dx;
    d2Cdx2(i) = (C(i+1)-2.*C(i)+C(i-1))./(dx.^2);
    dCdt(i) = D_S.*d2Cdx2(i) - k_E.*C(i);
end
for i = ns %right before interface - stroma
    dCdx(i) = (C_intfa - C(i))./dx;
    d2Cdx2(i) = (C(i-1)-2.*C(i)+C_intfa)./(dx.^2);
    dCdt(i) = D_S.*d2Cdx2(i) - k_E.*C(i);
end
for i = ns+1 %right after interface - OH
    dCdx(i) = (C(i+1)-C(i))./dx;
    d2Cdx2(i) = ((C_intfb)-2.*C(i)+C(i+1))./(dx.^2);
    dCdt(i) = D_OH.*d2Cdx2(i);
end
for i = (ns+2):(nOH+ns)-1 %in OH
    dCdx(i) = (C(i+1)-C(i))./dx;
    d2Cdx2(i) = (C(i+1)-2.*C(i)+C(i-1))./(dx.^2);
    dCdt(i) = D_OH.*d2Cdx2(i);
end  
for i = (nOH+ns)%end of stroma
    %dCdt(i) = D_S.*(-2*C(i)+C(i-1))./(dx.^2) - k_B.*C(i); %BC zero conc CORRECT
    dCdt(i) = 2*D_OH*(C(i-1)-C(i))./(dx.^2); %BC zero flux
end

DCDT = dCdt';

end

function S = calcsparsity(x,numx,nOH,ns)
    totalSize = numx;
    S = ones(totalSize,totalSize);
end