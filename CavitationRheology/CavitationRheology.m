t = 0:0.05:10*60; %time in seconds. 
Q_needle = 1/(3600); %1ml/hr converted to ml/sec;
E = 10^(4)*10; %elastic modulus of cervical tissue. In pascals, converted to dyn/cm2. From "Mechanical Properties of Female Reproductive Organs and Supporting Connective Tissues: A Review of the Current State of Knowledge"
mu = 1.0995*10^(-2); %cP viscosity ethanol, converted here to Poise. From "Density, viscosity and surface tension of water+ethanol mixtures from 293k to 323k".
L_needle = 31.75/10; %length of 27 gauge needle 1 1/4", or 31.75mm in length, changed here to cm. 
r_needle = 0.2032/(2*10); %nominal inner diameter of 27 gauge needle is 0.2032mm. Here we change to cm. 
r0 = r_needle; %inital radius of cavity = needle radius.

r = r0*(Q_needle*(t/((4/3)*pi*(r0^3)))+1).^(1/3);
lambda = r/r0;

Res_needle = (8*mu*L_needle)/(pi*(r_needle^4)); %resistance in needle
Pc = E*(5/6-(2./(3*lambda))-(1./(6*(lambda.^4)))); %pressure in cavity
Pup = Res_needle*Q_needle + Pc; %pressure upstream of needle
Pcrit = (5/6)*E*ones(size(lambda));
Stress = E*(lambda.^2-(1./lambda));

figure()
plot(t/60,r,'LineWidth',2,'color','k')
xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
ylabel('Cavity Radius (cm)','FontSize',18,'FontWeight','Bold')

figure()
plot(t/60,lambda,'LineWidth',2,'color','k')
xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
ylabel('\lambda','FontSize',18,'FontWeight','Bold')

figure()
plot(lambda,Pc,'LineWidth',2,'color','k')
hold on
plot(lambda,Pcrit,'LineWidth',2,'color','k','LineStyle','--')
%plot(lambda,Pup,'LineWidth',2,'color','r')

xlabel('\lambda','FontSize',18,'FontWeight','Bold')
ylabel('P (dyn/cm^{2})','FontSize',18,'FontWeight','Bold')
legend('P_{c}','P_{crit}')

figure()
plot(lambda,Pup,'LineWidth',2,'color','k')
xlabel('\lambda','FontSize',18,'FontWeight','Bold')
ylabel('P_{up} (dyn/cm^{2})','FontSize',18,'FontWeight','Bold')

figure()
plot(lambda,Stress,'LineWidth',2,'color','k')
xlabel('\lambda','FontSize',18,'FontWeight','Bold')
ylabel('\sigma (dyn/cm^{2})','FontSize',18,'FontWeight','Bold')
