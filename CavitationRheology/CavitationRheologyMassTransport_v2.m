clear
clc

%% DISCRETIZATION PARAMS

%t = 0:0.1:6*60; %time in seconds. 
t = linspace(0,5*60*60,1201);
Rtot = 2; %in cm, total dimensions of complete analysis

numr = 1001; %number of spatial steps
r = linspace(0,Rtot,numr); %radial distance from center
% dx = 0.005;
% r = 0:dx:Rtot; %radial distance from center
% numr = numel(r);

%% INJECTION PARAMS

Q_needle = 1/(3600); %1ml/hr converted to ml/sec;
L_needle = 31.75/10; %length of 27 gauge needle 1 1/4", or 31.75mm in length, changed here to cm. 
r_needle = 0.2032/(2*10); %nominal inner diameter of 27 gauge needle is 0.2032mm. Here we change to cm. 
Vol = 0.10; %total injection volume in ml
t_injection_end = Vol/Q_needle;
additional_drug = 0; % flag to determine if ethanol has been injected with an additional drug or not.
%r_needle = 0.01;

%% ETHANOL PARAMS

D_S = 7.5*10^(-6); %(cm^2/s)
D_C = D_S;
C_0 = 1; %normalized concentration of ethanol
phi = 1; %partition coefficient between cavity and tissue.
mu_oh = 1.0995*10^(-3)*1000/100; %cP viscosity ethanol, converted here to Poise. From "Density, viscosity and surface tension of water+ethanol mixtures from 293k to 323k".
mu_wat = 0.8914*10^(-3)*1000/100; %cP viscosity water, converted 
%% DRUG PARAMS

D_Sdrug = 5*10^(-7); %(cm^2/s)
D_Cdrug = D_Sdrug;
C_0drug = 1; %normalized concentration of ethanol
phidrug = 1; %partition coefficient between cavity and tissue.

%% HOST PARAMS

E = 10^(4)*10; %elastic modulus of cervical tissue. In pascals, converted to dyn/cm2. From "Mechanical Properties of Female Reproductive Organs and Supporting Connective Tissues: A Review of the Current State of Knowledge"
r0 = r_needle; %inital radius of cavity = needle radius.
K_wat = 1*10^(-14)*((100^4)/(10^5)); %hydraulic permeability in (m^4)/(N*s), here converted to cm^4/(dyne*s) = (cm^3 s^2)/g from "Direct Measurement of the Permeability of Human Cervical Tissue
K_oh = K_wat * mu_oh / mu_wat;
%k = K*(0.001/0.001095); %hydraulic conductivity in cervical tissue. In cm/s. CHANGE THIS. 
a0 = 1; %radius of tumor in cm.

%% SOLVER PARAMS

tspan = [0 t_injection_end];
n_initialcavityradius = find(r>=r0,1);
y0 = 1;
S = calcsparsity(r,numr);
S = sparse(S);
opts1 = odeset('Vectorized','on','JPattern',S);

%% VISUALIZATION PARAMS

show_strain_radius = 0;
show_pressures = 0;
show_heatmaps = 0;
show_concentration_lineplots = 1;
show_mass_conservation = 1;

%% CHANGING HYDRAULIC CONDUCTIVITY

Kvec = K_oh;
%tt = zeros(length(t),length(Kvec));
%lambda = zeros(length(t),length(Kvec));
%for i = 1:length(Kvec)
    
[tt_lambda,lambda] = ode45(@(tt,lambda) dLdt(tt,lambda,Q_needle,E,mu_oh,L_needle,r_needle,r0,K_oh,Rtot), tspan, y0);

nt_injection_end = find(t>=t_injection_end,1);

%tt = t_temp;
%lambda = lambda_temp;
%tt(:,i) = padarray(t_temp,size(tt,1)-length(t_temp),'post');
%lambda(:,i) = padarray(lambda_temp,size(lambda,1)-length(lambda_temp),'post');

%end

Res_needle = (8*mu_oh*L_needle)/(pi*(r_needle^4)); %resistance in needle
Pc = E*(5/6-(2./(3*lambda))-(1./(6*(lambda.^4)))); %pressure in cavity
Pup = Res_needle*Q_needle + Pc; %pressure upstream of needle
Pcrit = (5/6)*E*ones(size(lambda));
Stress = E*(lambda.^2-(1./lambda));

cavity_radius = lambda .* r0;

r_cavity = [tt_lambda,cavity_radius];


if (show_strain_radius == 1)
    plotCavityStrainRadius(tt_lambda, lambda, cavity_radius, 'k', Kvec)
end
if (show_pressures == 1)
    plotPressures(tt_lambda, lambda, Pc, Pup, Pcrit, 'k', Kvec)
end


nRetcell = zeros(1,length(t));

C_OH = zeros(length(tt_lambda),length(r));
C_drug = zeros(length(tt_lambda),length(r));

%ti = t(2:6);
%te = 0;

nRinit = find(r>=cavity_radius(1),1);
nRmax = find(r>=cavity_radius(end),1);
IC = zeros(1,numr); %Initial condition
IC(1:nRinit) = C_0;
IC_drug = zeros(1,numr);
IC_drug(1:nRinit) = C_0drug;

% v = zeros(1,length(r));
% dvdr = zeros(1,length(r));
%v = [0,Q_needle./(4*pi*r(2:end).^2)];
%dvdr = [0,-Q_needle./(2*pi*r(2:end).^3)];
% tmasstransport = zeros(1,length(ti));

%C_OH(1,:) = IC;
%C_drug(1,:) = IC;

tic 

% nr_tot = zeros(numel(tt),1);
% for i = 1:numel(tt)
%     [DCDT, nr_cavity] = dcdt(tt(i),0,r,D_S,D_C,phi,K_oh,a,Q_needle,r_cavity,Pc,numr);
%     nr_tot(i) = nr_cavity;
% end
% 
% figure(1)
% hold on
% plot(tt/60,r(nr_tot),'r')

[tt,C_OH] = ode15s(@(tt,C_OH) dcdt(tt,C_OH,r,D_S,D_C,phi,K_oh,Rtot,Q_needle,r_cavity,Pc,numr), tspan, IC, opts1);
IC_remain = C_OH(end,:);


if (t_injection_end < t(end))
    remainTspan = [tt(end)+tt(2)-tt(1), t(end)];
    
%     nr_tot = zeros(numel(tt),1);
%     for i = 1:numel(tt)
%         [DCDT, nr_cavity] = dcdt(remainT(i),Ct,r,D_S,D_C,phi,K_oh,Rtot,0,r_cavity,Pc,numr);
%         nr_tot(i) = nr_cavity;
%     end
% 
%     figure(1)
%     hold on
%     plot(tt/60,r(nr_tot),'r')
    
    [remainT,Ct] = ode15s(@(remainT,Ct) dcdt(remainT,Ct,r,D_S,D_C,phi,K_oh,Rtot,0,r_cavity,Pc,numr), remainTspan, IC_remain, opts1);
    C_OH = [C_OH; Ct];
    tt = [tt; remainT];
end

C_OH = C_OH';

toc

[nr_cavity, nt] = get_nr_cavity(r,r_cavity, tt);

% figure()
% plot(tt/60,r(nr_cavity),'LineWidth',2)



% for jj = 1:nt_injection_end
% 
%         %nRetcell = jj;
%         if jj == 1
% %             t_current = [0,ti];
%             %C(jj,:) = IC;
%             tspan = [tt(jj),tt(jj+2)];
%         elseif jj == nt_injection_end
%             tspan = [tt(jj-2),tt(jj)];
%         else
%             tspan = [tt(jj-1),tt(jj+1)];
%         end
%         
%         nRetcell(jj) = find(r>=cavity_radius(jj),1);
%         v = [0, Q_needle./(4*pi*r(2:nRetcell).^2), -K*a./(r(nRetcell+1:end).^2).*(Pc(jj)./(1-(a./r(nRetcell(jj)))))];
%         dvdr = [0,-Q_needle./(2*pi*r(2:end).^3), 2*K*a.*Pc(jj)./((r(nRetcell+1:end).^3).*(1-(a./r(nRetcell(jj)))))];
%                 
%         [t_current,Ct] = ode15s(@(t_current,Ct) dcdt(t_current,Ct,r,D_S,D_C,phi,numr,nRetcell(jj),v,dvdr,C_0), tspan, IC, opts1);
%         
% %         if (additional_drug == 1)
% %             [t_current_drug,Ct_drug] = ode15s(@(t_current_drug,Ct_drug) dcdt(t_current,Ct_drug,r,D_Sdrug,D_Cdrug,phidrug,numr,nRetcell(jj),v,dvdr,C_0drug), tspan, IC_drug, opts1);
% %             IC_drug = Ct_drug(end,:);
% % 
% %         end
%         IC = Ct(end,:);
%         if jj > 2 && jj < nt_injection_end
%             C_OH(jj,:) = mean(Ct,1);
% %             if (additional_drug == 1)
% %                 C_drug(jj,:) = mean(Ct_drug,1);
% %             end
%         elseif jj == nt_injection_end
%             C_OH(jj,:) = Ct(end,1);
%         end
%         
% 
% end

%% INJECTION END

% if (t_injection_end < t(end))
%     nRetcell(nt_injection_end:end) = nRetcell(nt_injection_end-1);
%     v = zeros(1,length(r));
%     dvdr = zeros(1,length(r));
%     remainTspan = [tt(end)+tt(2)-tt(1), t(end)];
%     [remainT,Ct] = ode15s(@(remainT,Ct) dcdt(remainT,Ct,r,D_S,D_C,phi,numr,nRetcell(nt_injection_end),v,dvdr,C_0), remainTspan, IC, opts1);
%     
% %     if (additional_drug == 1)
% %         [remainT_drug,Ct_drug] = ode15s(@(remainT_drug,Ct_drug) dcdt(remainT_drug,Ct_drug,r,D_Sdrug,D_Cdrug,phidrug,numr,nRetcell(nt_injection_end),v,dvdr,C_0drug), remainTspan, IC_drug, opts1);
% %         C_drug = [C_drug; Ct_drug];
% %     end
%     C_OH = [C_OH; Ct];
%     tt = [tt; remainT];
% % else 
% %     t_merged = tt;
% end

%%

% C_OH = C_OH';
% if (additional_drug == 1)
%     C_drug = C_drug';
% end

a = (a0^3-r0^3+r(nr_cavity).^3).^(1/3)';
a_vec = [tt,a];

tvec1=[find(tt>=0*60,1) find(tt>=1*60,1) find(tt>=2*60,1) find(tt>=4*60,1) find(tt>=6*60,1) find(tt>=8*60,1) find(tt>=12*60,1) find(tt>=20*60,1) find(tt>=30*60,1)];

%% LINE PLOTS ETHANOL VS R

if (show_concentration_lineplots == 1)
    plotConcentrationProfile(r, tvec1, C_OH)
end

%% SPHERE PLOTS 

if (show_heatmaps == 1)
    makeSphereHeatmaps(tt,tvec1,r,C_OH,a_vec)
end

%% LINE PLOTS MASS CONSERVATION VS TIME

if (show_mass_conservation == 1)
    %save test_massconservation_v2 r Q_needle Vol tt C_OH C_0 a
    plotNormalizedMassConservation(r, Q_needle, Vol, tt, C_OH, C_0, a)
end

%% approximation of cavity radius vs. deterministic result:

figure()
plot(tt_lambda/60,cavity_radius,'LineWidth',2)
hold on
plot(tt/60,r(nr_cavity),'k','LineWidth',2)
set(gca,'FontSize',18)
xlim([0,30])
xlabel('Time (mins)', 'FontSize', 22)
ylabel('Cavity Radius (cm)', 'FontSize', 22)
legend('computed','interpolated')


function [nr_cavity, nt] = get_nr_cavity(r,r_cavity, t)

    nr_cavity = zeros(length(t),1);
    nt = zeros(length(t),1);

    for i = 1:length(t)
        nt(i) = binarySearchBin(r_cavity(:,1),t(i));
        if (nt(i) == -1)
            nt(i) = numel(r_cavity(:,1));
        end
        nr_cavity(i) = binarySearchBin(r,r_cavity(nt(i),2));

    end

end

function [DLAMBDADT] = dLdt(t,lambda,Q_needle,E,mu,L_needle,r_needle,r0,k,a)

        Pc = E*(5/6-(2./(3*lambda))-(1./(6*(lambda.^4))));
        DLAMBDADT = (Q_needle./(4*pi*(r0^3)*(lambda.^2)))+((k./((r0^3).*(lambda.^2))).*(a./(1-(a./(lambda*r0))))*Pc);

end

%function DCDT = dcdt(t,C,r,D_S,D_C,phi,numr,nRetcell,v,dvdr,C_0)
function [DCDT] = dcdt(t,C,r,D_S,D_C,phi,K,a,Q,r_cavity,Pc,numr)

    %NOTE: TO UNDERSTAND THESE EQUATIONS, LOOK UP FOR DIFFERENCE EQUATIONS: CENTRAL, BACKWARD
    %AND FORWARD.

    dr = r(2)-r(1);
    dCdt = zeros(length(r),size(C,2));
    dCdr = zeros(length(r),size(C,2));
    d2Cdr2 = zeros(length(r),size(C,2));
    
    %         nRetcell(jj) = find(r>=cavity_radius(jj),1);
    
%     if (Q == 0)
%         nt = numel(r_cavity(:,1));
%         nr_cavity = nRmax;
%     else

%     nt = binarySearchBin(r_cavity(:,1),t);
%     if (nt == -1)
%         nt = numel(r_cavity(:,1));
%     end
%     nr_cavity = binarySearchBin(r,r_cavity(nt,2));

      [nr_cavity, nt] = get_nr_cavity(r,r_cavity, t);
    
    
%     end
    
%     DCDT = 0;
%     return
    
    %nr_cavity = find(r>=r_cavity(nt,2),1);

    %C at cavity/stroma interface
    C_intf_CSa = (D_C.*C(nr_cavity)+(D_S.*C(nr_cavity+1)))./((D_S.*phi)+D_C); %Right before interface
    C_intf_CSb = (D_C.*C(nr_cavity)+(D_S.*C(nr_cavity+1)))./((D_C./phi)+D_S); %Right after interface
    
    %v = [0,Q./(4*pi*r(2:end).^2)];
    %dvdr = [0,-Q./(2*pi*r(2:end).^3)];
   
    %v = [0, Q./(4*pi*r(2:nr_cavity).^2), -K*a./(r(nr_cavity+1:end).^2).*(Pc(nt)./(1-(a./r(nr_cavity))))];
    %dvdr = [0,-Q./(2*pi*r(2:end).^3), 2*K*a.*Pc(nt)./((r(nr_cavity+1:end).^3).*(1-(a./r(nr_cavity))))];
    
    %v = zeros(length(r));
    %dvdr = zeros(length(r));

    for i = 1 %y direction BC zero flux (dcdx=0)
        dCdt(i,:) = 6*D_C.*(C(i+1,:)-C(i,:))./(dr.^2);
    end
    for i = 2:nr_cavity-1 %inside cavity
        v = Q./(4*pi*r(i).^2);
        dvdr = -Q./(2*pi*r(i).^3);
        dCdr(i,:) = (C(i,:)-C(i-1,:))./dr; %backward difference
        %dCdr(i,:) = (C(i+1,:)-C(i-1,:))./(2*dr); %central difference
        %dCdr(i) = (C(i+1)-C(i))./dr; %forward difference
        d2Cdr2(i,:) = (C(i+1,:)-2.*C(i,:)+C(i-1,:))./(dr.^2);
        dCdt(i,:) = ((2*D_C/r(i)).*dCdr(i,:))+(D_C*d2Cdr2(i,:))-((2/r(i)).*(v.*C(i,:)))-(C(i,:).*dvdr)-(v.*dCdr(i,:));
        %dCdt(i,:) = ((2*D_C/r(i)).*dCdr(i,:))+(D_C*d2Cdr2(i,:))-((2/r(i)).*(v(i).*C(i,:)))-(C(i,:).*dvdr(i))-(v(i).*dCdr(i,:));
    end
    for i = nr_cavity %right before interface - cavity/stroma
        v = Q./(4*pi*r(i).^2);
        dvdr = -Q./(2*pi*r(i).^3);
        dCdr(i,:) = (C(i,:)-C(i-1,:))./dr; %backward difference
        %dCdr(i,:) = (C_intf_CSa-C(i-1,:))./(2*dr); %central difference
        %dCdr(i) = (C_intf_CSa-C(i))./dr; %forward difference
        d2Cdr2(i,:) = (C(i-1,:)-2.*C(i,:)+C_intf_CSa)./(dr.^2);
        dCdt(i,:) = ((2*D_C/r(i)).*dCdr(i,:))+(D_C*d2Cdr2(i,:))-((2/r(i)).*(v.*C(i,:)))-(C(i,:).*dvdr)-(v.*dCdr(i,:));
        %dCdt(i,:) = ((2*D_C/r(i)).*dCdr(i,:))+(D_C*d2Cdr2(i,:))-((2/r(i)).*(v(i).*C(i,:)))-(C(i,:).*dvdr(i))-(v(i).*dCdr(i,:));

    end
    for i = (nr_cavity+1) %right after interface - cavity/stroma
        %v = Q./(4*pi*r(i).^2);
        %dvdr = -Q./(2*pi*r(i).^3);
        v = (Q ~=0) .* -K*a./(r(i).^2).*(Pc(nt)./(1-(a./r(nr_cavity))));
        dvdr = (Q ~=0) .* 2*K*a.*Pc(nt)./((r(i).^3).*(1-(a./r(nr_cavity))));
        
        dCdr(i,:) = (C(i,:)-C_intf_CSb)./dr; %backward difference
        %dCdr(i,:) = (C(i+1,:)-C_intf_CSb)./(2*dr); %central difference
        %dCdr(i) = (C(i+1)-C(i))./dr; %forward difference
        d2Cdr2(i,:) = ((C_intf_CSb)-2.*C(i,:)+C(i+1,:))./(dr.^2);
        dCdt(i,:) = ((2*D_S/r(i)).*dCdr(i,:))+(D_S*d2Cdr2(i,:))-((2/r(i)).*(v.*C(i,:)))-(C(i,:).*dvdr)-(v.*dCdr(i,:));
        %dCdt(i,:) = ((2*D_S/r(i)).*dCdr(i,:))+(D_S*d2Cdr2(i,:))-((2/r(i)).*(v(i).*C(i,:)))-(C(i,:).*dvdr(i))-(v(i).*dCdr(i,:));

    end
    for i = nr_cavity+2:(numr-1) %in stroma
        %v = Q./(4*pi*r(i).^2);
        %dvdr = -Q./(2*pi*r(i).^3);
        v = (Q ~=0) .* -K*a./(r(i).^2).*(Pc(nt)./(1-(a./r(nr_cavity))));
        dvdr = (Q ~=0) .* 2*K*a.*Pc(nt)./((r(i).^3).*(1-(a./r(nr_cavity))));
        
        %dCdr(i) = (C(i)-C(i-1))./dr; %backward difference
        dCdr(i,:) = (C(i+1,:)-C(i-1,:))./(2*dr); %central difference
        %dCdr(i) = (C(i+1)-C(i))./dr; %forward difference
        d2Cdr2(i,:) = (C(i+1,:)-2.*C(i,:)+C(i-1,:))./(dr.^2);
        dCdt(i,:) = ((2*D_S/r(i)).*dCdr(i,:))+(D_S*d2Cdr2(i,:))-((2/r(i)).*(v.*C(i,:)))-(C(i,:).*dvdr)-(v.*dCdr(i,:));
        %dCdt(i,:) = ((2*D_S/r(i)).*dCdr(i,:))+(D_S*d2Cdr2(i,:))-((2/r(i)).*(v(i).*C(i,:)))-(C(i,:).*dvdr(i))-(v(i).*dCdr(i,:));

    end
    for i = numr %end of stroma, B.C.
        v = (Q ~=0) .* -K*a./(r(i).^2).*(Pc(nt)./(1-(a./r(nr_cavity))));
        
        dCdt(i,:) = D_S*(C(i,:)-2.*C(i-1,:)+C(i-2,:))./(dr.^2)-(2/a)*v*C(i);
        %dCdt(i,:) = D_S*(C(i,:)-2.*C(i-1,:)+C(i-2,:))./(dr.^2)-(2/a)*v(i)*C(i);

        %dCdt(i,:) = 6*D_S*(C(i-1,:)-C(i,:))./(dr.^2);
        %dCdt(i,:) = 2*D_S*(C(i-1,:)-C(i,:))./(dr.^2);
    end

    DCDT = dCdt;

end

 
function S = calcsparsity(r,numr)

    totalSize = numr;
    % vec = ones(3,totalSize);
    e = ones(totalSize,1);
    A = spdiags([e e e],-1:1,totalSize,totalSize);
    S = full(A);
    
    %S = ones(totalSize,totalSize);

end

