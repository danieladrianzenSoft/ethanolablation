
%% DISCRETIZATION PARAMS

t = 0:0.1:6*60; %time in seconds. 
Rtot = 2; %in cm, total dimensions of complete analysis
numr = 1601; %number of spatial steps
r = linspace(0,Rtot,numr); %radial distance from center

%% INJECTION PARAMS

Q_needle = 1/(3600); %1ml/hr converted to ml/sec;
L_needle = 31.75/10; %length of 27 gauge needle 1 1/4", or 31.75mm in length, changed here to cm. 
r_needle = 0.2032/(2*10); %nominal inner diameter of 27 gauge needle is 0.2032mm. Here we change to cm. 

%% DRUG PARAMS

D_S = 7.5*10^(-6); %(cm^2/s)
D_C = D_S;
C_0 = 1; %normalized concentration of ethanol
phi = 1; %partition coefficient between cavity and tissue.

%% HOST PARAMS

E = 10^(4)*10; %elastic modulus of cervical tissue. In pascals, converted to dyn/cm2. From "Mechanical Properties of Female Reproductive Organs and Supporting Connective Tissues: A Review of the Current State of Knowledge"
mu = 1.0995*10^(-2); %cP viscosity ethanol, converted here to Poise. From "Density, viscosity and surface tension of water+ethanol mixtures from 293k to 323k".
r0 = r_needle; %inital radius of cavity = needle radius.
K = 1*10^(-14)*((100^4)/(10^5)); %hydraulic permeability in (m^4)/(N*s), here converted to cm^4/(dyne*s) = (cm^3 s^2)/g from "Direct Measurement of the Permeability of Human Cervical Tissue
%k = K*(0.001/0.001095); %hydraulic conductivity in cervical tissue. In cm/s. CHANGE THIS. 
a = 1; %radius of tumor in cm.

%% SOLVER PARAMS
tspan = [0 t(end)];
n_initialcavityradius = find(r>=r0,1);
y0 = [1,C_0*ones(1,n_initialcavityradius),zeros(1,numr-n_initialcavityradius)];
S = calcsparsity(r,numr);
S = sparse(S);
opts1 = odeset('JPattern',S);


%% CHANGING HYDRAULIC CONDUCTIVITY

%Kvec = K;
%tt = zeros(length(t),1);
%lambda = zeros(length(t),length(Kvec));

%for i = 1:length(Kvec)
    
[tt,Y] = ode45(@(t,Y) dYdt(t,Y,r,Q_needle,E,mu,L_needle,r_needle,r0,K,a,D_C,D_S,phi,numr), t, y0, opts1);
%tt(:,i) = padarray(t_temp,size(tt,1)-length(t_temp),'post');
%lambda(:,i) = padarray(lambda_temp,size(lambda,1)-length(lambda_temp),'post');
%end
lambda = Y(:,1);
C = Y(:,2:end)';

Res_needle = (8*mu*L_needle)/(pi*(r_needle^4)); %resistance in needle
Pc = E*(5/6-(2./(3*lambda))-(1./(6*(lambda.^4)))); %pressure in cavity
Pup = Res_needle*Q_needle + Pc; %pressure upstream of needle
Pcrit = (5/6)*E*ones(size(lambda));
Stress = E*(lambda.^2-(1./lambda));


cavity_radius = lambda .* r0;

plotCavityStrainRadius(tt, lambda, cavity_radius, 'k', K)
plotPressures(tt, lambda, Pc, Pup, Pcrit, 'k', K)
ntumorbound = find(r>=a,1); %Tumor boundary at 0.5cm radius (1cm diameter)
tvec1=[find(tt>=0.5*60,1) find(tt>=1*60,1) find(tt>=3*60,1)];

makeSphereHeatmaps(tt,tvec1,r,C,ntumorbound)


function DYDT = dYdt(t,Y,r,Q_needle,E,mu,L_needle,r_needle,r0,k,a,D_C,D_S,phi,numr, tspan, y0)
        
        dr = r(2)-r(1);
        dYdt = zeros(length(r),size(Y,2));
        dYdr = zeros(length(r),size(Y,2));
        d2Ydr2 = zeros(length(r),size(Y,2));

        for i = 1
            Pc = E*(5/6-(2./(3*Y(i,:)))-(1./(6*(Y(i,:).^4))));
            dYdt(i,:) = (Q_needle/(4*pi*(r0^3)*(Y(i,:)^2)))+((k/((r0^3)*(Y(i,:)^2)))*(a/(1-(a/(Y(i,:)*r0))))*Pc);
        end
        
        nRetcell = find(r>=Y(1,:).*r0,1);
        
        %C at cavity/stroma interface
        C_intf_CSa = (D_C.*Y(nRetcell,:)+(D_S.*Y(nRetcell+1,:)))./((D_S.*phi)+D_C); %Right before interface
        C_intf_CSb = (D_C.*Y(nRetcell,:)+(D_S.*Y(nRetcell+1,:)))./((D_C./phi)+D_S); %Right after interface
        
        for i = 2 %y direction BC zero flux (dcdx=0)
            dYdt(i,:) = 6*D_C.*(Y(i+1,:)-Y(i,:))./(dr.^2);
        end
        for i = 3:nRetcell-1 %inside cavity
            dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr;
            d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
            dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
        end
        for i = nRetcell %right before interface - cavity/stroma
           dYdr(i,:) = (C_intf_CSa - Y(i,:))./dr;
           d2Ydr2(i,:) = (Y(i-1,:)-2.*Y(i,:)+C_intf_CSa)./(dr.^2);
           dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
        end
        for i = (nRetcell+1) %right after interface - cavity/stroma
            dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr;
            d2Ydr2(i,:) = ((C_intf_CSb)-2.*Y(i,:)+Y(i+1,:))./(dr.^2);
            dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:));
        end
        for i = nRetcell+2:(numr) %in stroma
            dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr;
            d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
            dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:));
        end
        for i = numr+1 %end of stroma, B.C.
            dYdt(i,:) = 6*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2);
        end

        DYDT = dYdt;

end

function S = calcsparsity(r,numr)

    totalSize = numr + 1;
    S = ones(totalSize,totalSize);

end
