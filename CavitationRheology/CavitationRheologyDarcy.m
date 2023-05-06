
t = 0:0.1:200; %time in seconds. 
Q_needle = 1/(3600); %1ml/hr converted to ml/sec;
E = 10^(4)*10; %elastic modulus of cervical tissue. In pascals, converted to dyn/cm2. From "Mechanical Properties of Female Reproductive Organs and Supporting Connective Tissues: A Review of the Current State of Knowledge"
mu = 1.0995*10^(-2); %cP viscosity ethanol, converted here to Poise. From "Density, viscosity and surface tension of water+ethanol mixtures from 293k to 323k".
L_needle = 31.75/10; %length of 27 gauge needle 1 1/4", or 31.75mm in length, changed here to cm. 
r_needle = 0.2032/(2*10); %nominal inner diameter of 27 gauge needle is 0.2032mm. Here we change to cm. 
r0 = r_needle; %inital radius of cavity = needle radius.
K = 1*10^(-14)*((100^4)/(10^5)); %hydraulic permeability in (m^4)/(N*s), here converted to cm^4/(dyne*s) = (cm^3 s^2)/g from "Direct Measurement of the Permeability of Human Cervical Tissue
%k = K*(0.001/0.001095); %hydraulic conductivity in cervical tissue. In cm/s. CHANGE THIS. 
a = 1; %radius of tumor in cm.

tspan = [0 t(end)];
y0 = 1;


%% CHANGING HYDRAULIC CONDUCTIVITY

Kvec = [K,K*10,K*100];
tt = zeros(length(t),length(Kvec));
lambda = zeros(length(t),length(Kvec));

for i = 1:length(Kvec)
    
[t_temp,lambda_temp] = ode45(@(t,lambda) dLdt(t,lambda,Q_needle,E,mu,L_needle,r_needle,r0,Kvec(i),a), t, y0);
tt(:,i) = t_temp;
lambda(:,i) = lambda_temp;
%tt(:,i) = padarray(t_temp,size(tt,1)-length(t_temp),'post');
%lambda(:,i) = padarray(lambda_temp,size(lambda,1)-length(lambda_temp),'post');

end

Res_needle = (8*mu*L_needle)/(pi*(r_needle^4)); %resistance in needle
Pc = E*(5/6-(2./(3*lambda))-(1./(6*(lambda.^4)))); %pressure in cavity
Pup = Res_needle*Q_needle + Pc; %pressure upstream of needle
Pcrit = (5/6)*E*ones(size(lambda));
Stress = E*(lambda.^2-(1./lambda));

fprintf('\nChanging k: \n')
fprintf('k: %d %d \n', size(Kvec))
fprintf('lambda: %d %d \n', size(lambda))
fprintf('Pcrit: %d %d \n', size(Pcrit))
fprintf('Pc: %d %d \n', size(Pc))

cavity_radius = lambda .* r0;

% arrayLegend = cell(1,length(k));
% for i = 1:length(k)
%     arrayLegend{i} = sprintf('k = %.0d cm^{4} dyn^{-1} s^{-1}',k(i));
% end

%cm(1,:) = [1,1,1]/255;
%cm(2,:) = [255,1,1]/255;
%cm(3,:) = [1,1,255]/255;

plotCavityStrainRadius(tt, lambda, cavity_radius, 'k', Kvec)
plotPressures(tt, lambda, Pc, Pup, Pcrit, 'k', Kvec)


%% CHANGING TISSUE ELASTICITY

E = 10^(4)*10;
Evec = [0.1*E, E, 10*E];

% tspan = [0 t(end)];
% y0 = 1;
tt = zeros(length(t),length(Evec));
lambda = zeros(length(t),length(Evec));

for i = 1:length(Evec)
    
[t_temp,lambda_temp] = ode45(@(t,lambda) dLdt(t,lambda,Q_needle,Evec(i),mu,L_needle,r_needle,r0,Kvec(1),a), t, y0);
tt(:,i) = t_temp;
lambda(:,i) = lambda_temp;
%tt(:,i) = padarray(t_temp,size(tt,1)-length(t_temp),'post');
%lambda(:,i) = padarray(lambda_temp,size(lambda,1)-length(lambda_temp),'post');

end

Res_needle = (8*mu*L_needle)/(pi*(r_needle^4)); %resistance in needle
Pc = (5/6-(2./(3*lambda))-(1./(6*(lambda.^4)))).*Evec; %pressure in cavity
Pup = Res_needle*Q_needle + Pc; %pressure upstream of needle
Pcrit = (5/6)*ones(size(lambda)).*Evec;
Stress = (lambda.^2-(1./lambda)).*Evec;

fprintf('\nChanging E: \n')
fprintf('Evec: %d %d \n', size(Evec))
fprintf('lambda: %d %d \n', size(lambda))
fprintf('Pcrit: %d %d \n', size(Pcrit))
fprintf('Pc: %d %d \n', size(Pc))

cavity_radius = lambda.*r0;
% 
% arrayLegend = cell(1,length(Evec));
% for i = 1:length(Evec)
%     arrayLegend{i} = sprintf('E = %.0e\t dyn cm^{-2}',Evec(i));
% end
% 
% figure()
% xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
% ylabel('\lambda','FontSize',18,'FontWeight','Bold')
% for i = 1:length(Evec)
% plot(tt(:,i)/60,lambda(:,i),'LineWidth',2,'Color',cm(i,:))
% hold on
% end
% legend(arrayLegend)
% set(gca,'FontSize',14)
% 
% figure()
% xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
% ylabel('Cavity Radius (cm)','FontSize',18,'FontWeight','Bold')
% for i = 1:length(Evec)
% plot(tt(:,i)/60,cavity_radius(:,i),'LineWidth',2,'Color',cm(i,:))
% hold on
% end
% legend(arrayLegend)
% set(gca,'FontSize',14)

plotCavityStrainRadius(tt, lambda, cavity_radius, 'E', Evec)
plotPressures(tt, lambda, Pc, Pup, Pcrit, 'E', Evec)

% figure()
% plot(lambda(:,3),Pc,'LineWidth',2,'color','k')
% hold on
% plot(lambda(:,3),Pcrit,'LineWidth',2,'color','k','LineStyle','--')
% plot(lambda(:,3),Pup,'LineWidth',2,'color','r')
%  
% xlabel('\lambda','FontSize',18,'FontWeight','Bold')
% ylabel('P (dyn/cm^{2})','FontSize',18,'FontWeight','Bold')
% legend('P_{c}','P_{crit}','P_{up}')
% set(gca,'FontSize',14)

%k=K*100;
%[t,lambda] = ode45(@(t,lambda) dLdt(t,lambda,Q_needle,E,mu,L_needle,r_needle,r0,k,a), tspan, y0);

% Res_needle = (8*mu*L_needle)/(pi*(r_needle^4)); %resistance in needle
% Pc = E*(5/6-(2./(3*lambda))-(1./(6*(lambda.^4)))); %pressure in cavity
% Pup = Res_needle*Q_needle + Pc; %pressure upstream of needle
% Pcrit = (5/6)*E*ones(size(lambda));
% Stress = E*(lambda.^2-(1./lambda));
% hold on
% plot(t/60,lambda,'LineWidth',2,'color','r')
% xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
% ylabel('\lambda','FontSize',18,'FontWeight','Bold')


% 
% figure() 
% plot(lambda,Pup,'LineWidth',2,'color','k')
% xlabel('\lambda','FontSize',18,'FontWeight','Bold')
% ylabel('P_{up} (dyn/cm^{2})','FontSize',18,'FontWeight','Bold')
% 
% figure()
% plot(lambda,Stress,'LineWidth',2,'color','k')
% xlabel('\lambda','FontSize',18,'FontWeight','Bold')
% ylabel('\sigma (dyn/cm^{2})','FontSize',18,'FontWeight','Bold')


function DLAMBDADT = dLdt(t,lambda,Q_needle,E,mu,L_needle,r_needle,r0,k,a, tspan, y0)

        Pc = E*(5/6-(2./(3*lambda))-(1./(6*(lambda.^4))));
        DLAMBDADT = (Q_needle/(4*pi*(r0^3)*(lambda^2)))+((k/((r0^3)*(lambda^2)))*(a/(1-(a/(lambda*r0))))*Pc);

end
