clear
clc

%% INJECTION PARAMS

Q_needle = 10/(3600); %1ml/hr converted to ml/sec; -> CHANGE TO 10 ML/HR for humans
L_needle = 31.75/10; %length of 27 gauge needle 1 1/4", or 31.75mm in length, changed here to cm. 
%r_needle = 0.21/(2*10); %nominal inner diameter of 27 gauge needle is 0.2032mm. Here we change to cm. 
r_needle = round(0.311/(2*10),3); %nominal inner diameter of 24 gauge needle is 0.2032mm. Here we change to cm. 
Vol = 0.50; %total injection volume -> CHANGE TO 0.5ML - 2.5ML for humans
%Vol = 1;
%t_injection_end = Vol/Q_needle; % BASED ON HUMAN MODELS, WOULD VARY FROM 3 TO 15 MINUTES
additional_drug = 1; % flag to determine if ethanol has been injected with an additional drug or not.
omega = 0.03;   %ethyl cellulose concentration relative to EC
rho_etoh = 0.789; % density of ethanol.

%% DISCRETIZATION PARAMS

dt_injection = 1; %1 second
dt_relaxation = 6*60; %6 mins
t_injection_start = 0;
t_injection_end = Vol/Q_needle;
t_end = 1*12*60*60;

%t_injection = 0:dt_injection:t_injection_end;
%t_relaxation = t_injection_end+dt_injection:dt_relaxation:t_end;
%t = [t_injection, t_relaxation];
Rtot = 6; %in cm, total dimensions of complete analysis
%numr = 1001; %number of spatial steps
dx_r0 = 0.002;
dx_rend = 0.005;
%dx_rend = 0.0025;

%dx_r0 = 0.001;
%dx_r1 = 0.003;
%dx_r2 = 0.005;
%r = 0:dx:Rtot;
r = [0:dx_r0:r_needle, r_needle+dx_r0:dx_rend:Rtot]; %radial distance from center
numr = numel(r);

%% CHANGING OMEGA (ETHANOL CONCENTRATION)

phiVec = [1; 1-omega; 1-2*omega; 1-3*omega; 1-4*omega];
%phiVec = 1-2*omega;
%phiVec = [1; 1-omega; 1-2*omega; 1-4*omega];

%% ETHANOL PARAMS

D_S = 5*10^(-6); %(cm^2/s)
%D_S = 1*10^(-7); %(cm^2/s) TEST
%D_C = D_S*0.6;
C_0 = 1; %normalized concentration of ethanol
P = 1; %partition coefficient between cavity and tissue.
%phi_vec = linspace(0.6,1,length(omegaVec)); %assuming partition coefficient varies with ECE concentration. Higher ECE leads to higher partition.
mu_oh = 1.0995*10^(-3)*1000/100; %cP viscosity ethanol, converted here to Poise. From "Density, viscosity and surface tension of water+ethanol mixtures from 293k to 323k".
mu_wat = 0.8914*10^(-3)*1000/100; %cP viscosity water, converted 

%% DRUG PARAMS

D_Sdrug = 5*10^(-7); %(cm^2/s)
D_Cdrug = D_Sdrug*1.5;
C_0drug = 1; %normalized concentration of ethanol
Pdrug = 1; %partition coefficient between cavity and tissue.


%% HOST PARAMS (+ HYDRAULIC CONDUCTIVITY, DEPENDENT ON DRUG)

E = 10^(4)*10; %elastic modulus of cervical tissue. In pascals, converted to dyn/cm2. From "Mechanical Properties of Female Reproductive Organs and Supporting Connective Tissues: A Review of the Current State of Knowledge"
%r0 = r_needle; %inital radius of cavity = needle radius.
K_wat_t = 1*10^(-14)*((100^4)/(10^5)); %hydraulic permeability in (m^4)/(N*s), here converted to cm^4/(dyne*s) = (cm^3 s)/g from "Direct Measurement of the Permeability of Human Cervical Tissue
K_oh_t = K_wat_t * mu_oh / mu_wat;
K_drug_t = K_wat_t * mu_oh / mu_wat;
K_oh_c = 3*10^(-8)/mu_oh; % CHECK ON THIS Effect of ethanol on water permeability of controlled release films composed of ethyl
%k = K*(0.001/0.001095); %hydraulic conductivity in cervical tissue. In cm/s. CHANGE THIS. 
a0 = 1; %radius of tumor in cm.
epsilon_t = 0.3; % porosity in tissue
%kappa = 0.3; % percentage of water 'imbibed' by the gel.
kappa = 0;

%% CHANGING D_S, EPSILON_C (POROSITY IN CAVITY), P AND THETA DUE TO OMEGA

% vectors for D_C, epsilon_c and theta given 0, 3, 6, 9 and 12% EC concs.
D_C_allPhi = linspace(1*10^(-5), 5e-7, 5);
epsilon_c_allPhi = linspace(1, 0.8, 5);
theta_allPhi = epsilon_c_allPhi ./ epsilon_t;
%theta_allPhi = ones(length(phiVec),5);
P_allPhi = P.*ones(length(phiVec),5); % assuming partition coefficient is the same regardless of ECE concentration, and is = 1.

D_C_vec = zeros(numel(phiVec),1);
epsilon_c_vec = zeros(numel(phiVec),1);
theta_vec = zeros(numel(phiVec),1);
P_vec = zeros(numel(phiVec),1);

for i = 1:length(phiVec)
    switch phiVec(i)
        case 1
            D_C_vec(i) = D_C_allPhi(1);
            epsilon_c_vec(i) = epsilon_c_allPhi(1); 
            theta_vec(i) = theta_allPhi(1);
            P_vec(i) = P_allPhi(1);
        case 1 - omega
            D_C_vec(i) = D_C_allPhi(2);
            epsilon_c_vec(i) = epsilon_c_allPhi(2); 
            theta_vec(i) = theta_allPhi(2);
            P_vec(i) = P_allPhi(2);
        case 1 - 2*omega
            D_C_vec(i) = D_C_allPhi(3);
            epsilon_c_vec(i) = epsilon_c_allPhi(3); 
            theta_vec(i) = theta_allPhi(3);
            P_vec(i) = P_allPhi(3);
        case 1 - 3*omega
            D_C_vec(i) = D_C_allPhi(4);
            epsilon_c_vec(i) = epsilon_c_allPhi(4); 
            theta_vec(i) = theta_allPhi(4);
            P_vec(i) = P_allPhi(4);
        case 1 - 4*omega
            D_C_vec(i) = D_C_allPhi(5);
            epsilon_c_vec(i) = epsilon_c_allPhi(5); 
            theta_vec(i) = theta_allPhi(5);
            P_vec(i) = P_allPhi(5);
        otherwise
            fprintf('\nEC concentration not in considered values \n')
    end
end

%% SOLVER PARAMS

tspan = [0 Vol/Q_needle];
%n_initialcavityradius = find(r>=r0,1);
y0 = 1;
getSparsity = getSparsityLib;
%S = sparse(getSparsity.diagonalBandLastRowLastCol(numr));
S = sparse(getSparsity.diagonalBandLastElem(numr));
S = sparse(S);
%'MaxStep', 1e-1
%opts1 = odeset('Vectorized','on', 'JPattern', S, 'Stats', 'off', 'RelTol', 1e-11, 'AbsTol', 1e-12, 'MaxStep', 1e-1,'BDF','on','MaxOrder',4);
opts1 = odeset('Vectorized','on','JPattern',S,'Refine',1);

keySet = ["Q_needle","L_needle","r_needle", "Vol", "rho_etoh", ...
    "Rtot","D_S","numr","a0","epsilon_t", ...
    "kappa", "y0"];
valueSet = [Q_needle,  L_needle, r_needle,  Vol, rho_etoh, ...
           Rtot, D_S, numr,  a0,  epsilon_t, ...
           kappa, y0];

params = dictionary(keySet, valueSet);
options = dictionary('tolerance',1e-3);


%% VISUALIZATION PARAMS

show_strain = 0;
show_cavity_radius = 1;
show_pressure_field = 0;
show_velocity_field = 1;
%show_pc_pcrit_p0 = 0;
show_pc_pcrit = 1;
show_p_inf = 0;
show_heatmaps = 0;
show_concentration_lineplots = 1;
show_ethanol_cloud_volume = 1;
show_mass_conservation = 1;
show_cavity_volume = 1;

