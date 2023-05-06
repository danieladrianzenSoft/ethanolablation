t = 0:0.1:10*60; %time in seconds. 
Q_needle = 1/(3600); %1ml/hr converted to ml/sec;
E = 10^(4)*10; %elastic modulus of cervical tissue. In pascals, converted to dyn/cm2. From "Mechanical Properties of Female Reproductive Organs and Supporting Connective Tissues: A Review of the Current State of Knowledge"
mu = 1.0995*10^(-2); %cP viscosity ethanol, converted here to Poise. From "Density, viscosity and surface tension of water+ethanol mixtures from 293k to 323k".
L_needle = 31.75/10; %length of 27 gauge needle 1 1/4", or 31.75mm in length, changed here to cm. 
r_needle = 0.2032/(2*10); %nominal inner diameter of 27 gauge needle is 0.2032mm. Here we change to cm. 
r0 = r_needle; %inital radius of cavity = needle radius.
K = 1*10^(-14)*((100^4)/(10^5)); %hydraulic permeability in (m^4)/(N*s), here converted to cm^4/(N*s) =  from "Direct Measurement of the Permeability of Human Cervical Tissue
k = 1; %hydraulic conductivity in cervical tissue. In cm/s. CHANGE THIS. 
a = 1; %radius of tumor in cm.

lambda = ones(1,length(t));
Pc = ones(1,length(t));
%lambdaTemp = 1;
% Pc = E*(5/6-(2./(3*lambdaTemp))-(1./(6*(lambdaTemp.^4))));
% trapz(a*Pc/(1-(a/(lambdaTemp*r0))));

for i = 2
    diff = 1;
    t_temp = t(1:i);
    lambdaTemp = lambda(i);
    while diff>=0.1
        LHS = Q_needle*t(i);
        Pc_temp = E*(5/6-(2./(3*lambdaTemp))-(1./(6*(lambdaTemp.^4))));
        RHS = (4/3)*pi*(r0^(3))*((lambdaTemp.^3)-1)+4*pi*k*a*trapz(Pc_temp./(1-(a./(lambdaTemp*r0))),1);
        diff = abs(LHS-RHS)/(LHS)
        lambdaTemp = lambdaTemp + 0.1;
    end
    lambda(i) = lambdaTemp - 0.1;
end

figure()
plot(t/60,lambda,'LineWidth',2,'color','k')
xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
ylabel('\lambda','FontSize',18,'FontWeight','Bold')
%r = r0*(Q_needle*(t/((4/3)*pi*(r0^3)))+1).^(1/3);
%lambda = r/r0;

%Res_needle = (8*mu*L_needle)/(pi*(r_needle^4)); %resistance in needle
%Pc = E*(5/6-(2./(3*lambda))-(1./(6*(lambda.^4)))); %pressure in cavity