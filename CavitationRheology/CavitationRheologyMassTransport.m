clear
clc

%% DISCRETIZATION PARAMS

t = linspace(0,10*60*60,1201);
Rtot = 2; %in cm, total dimensions of complete analysis
%numr = 1001; %number of spatial steps
dx = 0.001;
r = 0:dx:Rtot; %radial distance from center
numr = numel(r);

%% INJECTION PARAMS

Q_needle = 1/(3600); %1ml/hr converted to ml/sec; -> CHANGE TO 10 ML/HR for humans
L_needle = 31.75/10; %length of 27 gauge needle 1 1/4", or 31.75mm in length, changed here to cm. 
r_needle = 0.2032/(2*10); %nominal inner diameter of 27 gauge needle is 0.2032mm. Here we change to cm. 
Vol = 0.10; %total injection volume -> CHANGE TO 0.5ML - 2.5ML for humans
t_injection_end = Vol/Q_needle; % BASED ON HUMAN MODELS, WOULD VARY FROM 3 TO 15 MINUTES
additional_drug = 0; % flag to determine if ethanol has been injected with an additional drug or not.

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

%% HOST PARAMS (+ HYDRAULIC CONDUCTIVITY, DEPENDENT ON DRUG)

E = 10^(4)*10; %elastic modulus of cervical tissue. In pascals, converted to dyn/cm2. From "Mechanical Properties of Female Reproductive Organs and Supporting Connective Tissues: A Review of the Current State of Knowledge"
r0 = r_needle; %inital radius of cavity = needle radius.
K_wat = 1*10^(-14)*((100^4)/(10^5)); %hydraulic permeability in (m^4)/(N*s), here converted to cm^4/(dyne*s) = (cm^3 s^2)/g from "Direct Measurement of the Permeability of Human Cervical Tissue
K_oh = K_wat * mu_oh / mu_wat;
K_drug = K_wat * mu_oh / mu_wat;
%k = K*(0.001/0.001095); %hydraulic conductivity in cervical tissue. In cm/s. CHANGE THIS. 
a0 = 1; %radius of tumor in cm.

%% SOLVER PARAMS

tspan = [0 t_injection_end];
%n_initialcavityradius = find(r>=r0,1);
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
% show_halfhalf_heatmap = 1;



%% CHANGING HYDRAULIC CONDUCTIVITY

Kvec = K_oh;

%% FINDING STRAIN RATE, PRESSURES AND CAVITY RADIUS
    
[tt,lambda] = ode45(@(tt,lambda) dLdt(tt,lambda,Q_needle,E,mu_oh,L_needle,r_needle,r0,K_oh,Rtot), tspan, y0);

nt_injection_end = find(tt>=t_injection_end,1);

Res_needle = (8*mu_oh*L_needle)/(pi*(r_needle^4)); %resistance in needle
Pc = E*(5/6-(2./(3*lambda))-(1./(6*(lambda.^4)))); %pressure in cavity
Pup = Res_needle*Q_needle + Pc; %pressure upstream of needle
Pcrit = (5/6)*E*ones(size(lambda));
Stress = E*(lambda.^2-(1./lambda));

cavity_radius = lambda .* r0;

%% PLOTTING STRAIN AND PRESSURE

if (show_strain_radius == 1)
    plotCavityStrainRadius(tt, lambda, cavity_radius, 'k', Kvec)
end
if (show_pressures == 1)
    plotPressures(tt, lambda, Pc, Pup, Pcrit, 'k', Kvec)
end
%% MASS TRANSPORT MODEL

C_OH = zeros(length(tt),length(r));
C_drug = zeros(length(tt),length(r));

nRinit = find(r>=cavity_radius(1),1);

nRetcell = zeros(1,length(tt)); % vector that will contain index in r of cavity radius
nRetcell(1) = nRinit;

IC_oh = zeros(1,numr); %Initial condition
IC_oh(1:nRinit) = C_0;
IC_drug = zeros(1,numr);
IC_drug(1:nRinit) = C_0drug;


C_OH(1,:) = IC_oh;
C_drug(1,:) = IC_drug;

tic 

%%% LOOPING THROUGH TIME VECTOR, AND SOLVING MASS TRANSPORT
%%% AT EACH TIME POINT DURING INJECTION

for jj = 2:nt_injection_end

%         if jj == 1
%             tspan = [tt(jj),tt(jj+2)];
        if jj == nt_injection_end
            tspan = [tt(jj-2),tt(jj)];
        else
            tspan = [tt(jj-1),tt(jj+1)];
        end
       
        
        nRetcell(jj) = binarySearchBin(r,cavity_radius(jj));
        
        v = [0, Q_needle./(4*pi*r(2:nRetcell(jj)).^2), -K_oh*Rtot./(r(nRetcell(jj)+1:end).^2).*(Pc(jj)./(1-(Rtot./r(nRetcell(jj)))))];
        dvdr = [0, -Q_needle./(2*pi*r(2:nRetcell(jj)).^3), 2*K_oh*Rtot.*Pc(jj)./((r(nRetcell(jj)+1:end).^3).*(1-(Rtot./r(nRetcell(jj)))))];        
        
        [t_current_oh,Ct_oh] = ode15s(@(t_current,Ct) dcdt(t_current,Ct,r,D_S,D_C,phi,a0,numr,nRetcell(jj),v,dvdr,C_0), tspan, IC_oh, opts1);
        
        IC_oh = Ct_oh(end,:);
        if jj >= 2 && jj < nt_injection_end
            C_OH(jj,:) = mean(Ct_oh,1);
        elseif (jj == nt_injection_end)
            C_OH(jj,:) = Ct_oh(end,:);
        end
        
        if (additional_drug == 1)
            v = [0, Q_needle./(4*pi*r(2:nRetcell(jj)).^2), -K_drug*Rtot./(r(nRetcell(jj)+1:end).^2).*(Pc(jj)./(1-(Rtot./r(nRetcell(jj)))))];
            dvdr = [0, -Q_needle./(2*pi*r(2:nRetcell(jj)).^3), 2*K_drug*Rtot.*Pc(jj)./((r(nRetcell(jj)+1:end).^3).*(1-(Rtot./r(nRetcell(jj)))))];  
            
            [t_current_drug,Ct_drug] = ode15s(@(t_current_drug,Ct_drug) dcdt(t_current_drug,Ct_drug,r,D_Sdrug,D_Cdrug,phidrug,a0,numr,nRetcell(jj),v,dvdr,C_0drug), tspan, IC_drug, opts1);
        
            IC_drug = Ct_drug(end,:);
            if jj >= 2 && jj < nt_injection_end
                C_drug(jj,:) = mean(Ct_drug,1);
            elseif (jj == nt_injection_end)
                C_drug(jj,:) = Ct_drug(end,:);
            end

        end
        
end

%%% SOLVING MASS TRANSPORT FOR REMAINING TIME AFTER 
%%% INJECTION END

if (t_injection_end < t(end))
    v = zeros(1,length(r));
    dvdr = zeros(1,length(r));
    remainTspan = [tt(end)+tt(2)-tt(1), t(end)];
    [remainT,Ct_oh] = ode15s(@(remainT,Ct) dcdt(remainT,Ct,r,D_S,D_C,phi,a0,numr,nRetcell(nt_injection_end-1),v,dvdr,C_0), remainTspan, IC_oh, opts1);
    nr_cavity = [nRetcell,nRetcell(nt_injection_end) * ones(1,numel(remainT))];

    C_OH = [C_OH; Ct_oh];
    ttot = [tt; remainT];
    
    if (additional_drug == 1)
        [remainT_drug,Ct_drug] = ode15s(@(remainT_drug,Ct_drug) dcdt(remainT_drug,Ct_drug,r,D_Sdrug,D_Cdrug,phidrug,a0,numr,nRetcell(nt_injection_end-1),v,dvdr,C_0drug), remainTspan, IC_drug, opts1);
        nr_cavity_drug = [nRetcell,nRetcell(nt_injection_end) * ones(1,numel(remainT_drug))];

        C_drug = [C_drug; Ct_drug];
        ttot_drug = [tt; remainT_drug];
    end


end

C_OH = C_OH';
C_drug = C_drug';


toc

%% GETTING TUMOR RADIUS OVER TIME

a = (a0^3-r0^3+r(nr_cavity).^3).^(1/3)';

a_vec = [ttot,a];

figure()
plot(ttot/60,a,'LineWidth',2)

%ntumorbound = find(r>=a0,1); %Tumor boundary at 0.5cm radius (1cm diameter)
tvec_OH=[find(ttot>=0*60,1) find(ttot>=1*60,1) find(ttot>=2*60,1) find(ttot>=4*60,1) find(ttot>=6*60,1) find(ttot>=8*60,1) find(ttot>=12*60,1) find(ttot>=20*60,1) find(ttot>=30*60,1)];

if (additional_drug == 1)
    a_drug = (a0^3-r0^3+r(nr_cavity_drug).^3).^(1/3)';
    a_vec_drug = [ttot_drug,a_drug];
    tvec_drug=[find(ttot_drug>=0*60,1) find(ttot_drug>=1*60,1) find(ttot_drug>=2*60,1) find(ttot_drug>=4*60,1) find(ttot_drug>=6*60,1) find(ttot_drug>=8*60,1) find(ttot_drug>=12*60,1) find(ttot_drug>=20*60,1) find(ttot_drug>=30*60,1)];

end

%% LINE PLOTS ETHANOL VS R

if (show_concentration_lineplots == 1)

    plotConcentrationProfile(r, tvec_OH, C_OH)
    if (additional_drug == 1)
        plotConcentrationProfile(r, tvec_drug, C_drug)
    end

end

%% LINE PLOTS MASS CONSERVATION VS TIME

if (show_mass_conservation)
    save test_massconservation r nr_cavity Q_needle Vol ttot C_OH C_0 a
    plotNormalizedMassConservation(r, Q_needle, Vol, ttot, C_OH, C_0, a)
end

%% SPHERE PLOTS 

if (show_heatmaps == 1)

    if (additional_drug == 1)
        
        makeHalfAndHalfSphereHeatmaps(ttot,tvec_OH,tvec_drug,r,C_OH,C_drug,a_vec) 
        %makeSphereHeatmaps(ttot_drug,tvec_drug,r,C_drug,a_vec_drug)
    else
        makeSphereHeatmaps(ttot,tvec_OH,r,C_OH,a_vec)
    end

end
% if (show_halfhalf_heatmap == 1 && additional_drug == 1)
%     
%     %save test_sphereplot ttot tvec_OH tvec_drug r C_OH C_drug a_vec
%     makeHalfAndHalfSphereHeatmaps(ttot,tvec_OH,tvec_drug,r,C_OH,C_drug,a_vec) 
% 
% end
%% approximation of cavity radius vs. computed result:

figure()
plot(tt/60,cavity_radius,'LineWidth',2)
hold on
plot(ttot/60,r(nr_cavity),'k','LineWidth',2)
set(gca,'FontSize',18)
xlim([0,30])
xlabel('Time (mins)', 'FontSize', 22)
ylabel('Cavity Radius (cm)', 'FontSize', 22)
legend('computed','interpolated')


function DLAMBDADT = dLdt(t,lambda,Q_needle,E,mu,L_needle,r_needle,r0,k,a, tspan, y0)

        Pc = E*(5/6-(2./(3*lambda))-(1./(6*(lambda.^4))));
        DLAMBDADT = (Q_needle./(4*pi*(r0^3)*(lambda.^2)))+((k./((r0^3).*(lambda.^2))).*(a./(1-(a./(lambda*r0))))*Pc);

end


function DCDT = dcdt(t,C,r,D_S,D_C,phi,a,numr,nRetcell,v,dvdr,C_0)

    %NOTE: TO UNDERSTAND THESE EQUATIONS, LOOK UP FOR DIFFERENCE EQUATIONS: CENTRAL, BACKWARD
    %AND FORWARD.

    dr = r(2)-r(1);
    dCdt = zeros(length(r),size(C,2));
    dCdr = zeros(length(r),size(C,2));
    d2Cdr2 = zeros(length(r),size(C,2));

    %C at cavity/stroma interface
    C_intf_CSa = (D_C.*C(nRetcell)+(D_S.*C(nRetcell+1)))./((D_S.*phi)+D_C); %Right before interface
    C_intf_CSb = (D_C.*C(nRetcell)+(D_S.*C(nRetcell+1)))./((D_C./phi)+D_S); %Right after interface

    for i = 1 %y direction BC zero flux (dcdx=0)
        dCdt(i,:) = 6*D_C.*(C(i+1,:)-C(i,:))./(dr.^2);
    end
    for i = 2:nRetcell-1 %inside cavity
        dCdr(i,:) = (C(i,:)-C(i-1,:))./dr; %backward difference
        %dCdr(i,:) = (C(i+1,:)-C(i-1,:))./(2*dr); %central difference
        %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward difference
        d2Cdr2(i,:) = (C(i+1,:)-2.*C(i,:)+C(i-1,:))./(dr.^2);
        dCdt(i,:) = ((2*D_C/r(i)).*dCdr(i,:))+(D_C*d2Cdr2(i,:))-((2/r(i)).*(v(i).*C(i,:)))-(C(i,:).*dvdr(i))-(v(i).*dCdr(i,:));
    end
    for i = nRetcell %right before interface - cavity/stroma
       dCdr(i,:) = (C(i,:)-C(i-1,:))./dr; %backward difference
       %dCdr(i,:) = (C_intf_CSa-C(i-1,:))./(2*dr); %central difference
       %dCdr(i,:) = (C_intf_CSa-C(i,:))./dr; %forward difference
       d2Cdr2(i,:) = (C(i-1,:)-2.*C(i,:)+C_intf_CSa)./(dr.^2);
       dCdt(i,:) = ((2*D_C/r(i)).*dCdr(i,:))+(D_C*d2Cdr2(i,:))-((2/r(i)).*(v(i).*C(i,:)))-(C(i,:).*dvdr(i))-(v(i).*dCdr(i,:));
    end
    for i = (nRetcell+1) %right after interface - cavity/stroma
        dCdr(i,:) = (C(i,:)-C_intf_CSb)./dr; %backward difference
        %dCdr(i,:) = (C(i+1,:)-C_intf_CSb)./(2*dr); %central difference
        %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward difference
        d2Cdr2(i,:) = ((C_intf_CSb)-2.*C(i,:)+C(i+1,:))./(dr.^2);
        dCdt(i,:) = ((2*D_S/r(i)).*dCdr(i,:))+(D_S*d2Cdr2(i,:))-((2/r(i)).*(v(i).*C(i,:)))-(C(i,:).*dvdr(i))-(v(i).*dCdr(i,:));
    end
    for i = nRetcell+2:(numr-1) %in stroma
        %dCdr(i) = (C(i)-C(i-1))./dr; %backward difference
        dCdr(i,:) = (C(i+1,:)-C(i-1,:))./(2*dr); %central difference
        %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward difference
        d2Cdr2(i,:) = (C(i+1,:)-2.*C(i,:)+C(i-1,:))./(dr.^2);
        dCdt(i,:) = ((2*D_S/r(i)).*dCdr(i,:))+(D_S*d2Cdr2(i,:))-((2/r(i)).*(v(i).*C(i,:)))-(C(i,:).*dvdr(i))-(v(i).*dCdr(i,:));
    end
    for i = numr %end of stroma, B.C.
        %dCdt(i,:) = 6*D_S*(C(i-1,:)-C(i,:))./(dr.^2);
        dCdt(i,:) = D_S*(C(i,:)-2.*C(i-1,:)+C(i-2,:))./(dr.^2)-(2/a)*v(i)*C(i);

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

