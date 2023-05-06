clear
clc

%% initialize params

variableParam = "K";
%K = [3.8e-11, 7.6e-11, 3.8e-10, 7.6e-10]; % Netti 2003, 7.6 m^2 / kPa s = 7.6 cm^2 / barye s = 7.6 cm^3 s / g
K = 7.6e-11;
%K = [8e-10, 8e-11, 8e-12]; %hydraulic permeability 1e-11 from "Direct Measurement of the Permeability of Human Cervical Tissue
r0 = 0.035; % Netti 2003, 0.035cm
%r0 = round(0.311/(2*10),3); %nominal inner diameter of 24 gauge needle. Here we change to cm.
Rtot = 1; % Netti 2003, cm
%Rtot = 4; % cm
%Q = 1e-4 / 60; % Netti 2003, 0.1 mm^3 / min = 1e-4 cm^3 / min converted to s.
%Q = 10/(3600); % 10ml/hr converted to ml/sec;
Q = (0.1/10^3)/60; % 0.1uL/min Netti 2003
lambda = 13.16 * 1e4; % Netti 2003, 13.16kPa = 13.16*1e4 Barye
mu = 6.58 * 1e4; % Netti 2003, 6.58kPa = 6.58*1e4 Barye
c = K * (lambda + 2 * mu);
%c = [52,26,13,0] * 1e4; %figure 6 netti 2003, kPa converted to Barye
phi = 0.2; % Netti 2003, Porosity, dimensionless
r = linspace(0.01,Rtot,1000)';
t = 0:1:6*60;
nmax = 1e2; % accuracy of infinite series solution

% f = figure;
% ax = axes('Parent',f);
% plot(ax, r,besseljs(0,r),'LineWidth',2);
% hold on
% plot(ax, r,besseljs(1,r),'LineWidth',2);
% plot(ax, r,besseljs(2,r),'LineWidth',2);

%% make parameter dictionary

paramNames = ["K","r0","Rtot","Q","lambda","mu","phi","nmax"];
paramIndices = 1:length(paramNames);
paramValues = {K, r0, Rtot, Q, lambda, mu, phi, nmax};
params = dictionary(paramNames, paramValues);
paramLabels = dictionary(paramIndices,paramNames);

%% calculate dilatation, solid displacement, cavity radius and pressure

variableParamArray = params(variableParam);
lengthVariableParam = length(variableParamArray{1});

e = cell(lengthVariableParam,1);
u = cell(lengthVariableParam,1);
rc = cell(lengthVariableParam,1);
p = cell(lengthVariableParam,1);
e_laplace = cell(lengthVariableParam,1);

for iter = 1:lengthVariableParam
    %var_iter = var(iter);
    paramsIterValues = getParamValues(params,paramLabels,variableParam,iter);
    paramsIter = dictionary(paramNames,paramsIterValues);
    tic
    e{iter} = getDilatation(r,t,paramsIter);
    u{iter} = getDisplacement(r,t,paramsIter);
    toc
    rc{iter} = r0 + u{iter}(1,:);
    p{iter} = e{iter} * (lambda + 2*mu);
    %laplace transform approach
    fun = @(s) getDilatationLaplace(r,s,paramsIter);
    %e_laplace{iter} = talbot_inversion(fun,t);
end

%% plot results

% plot cavity radius vs time
f = figure;
ax = axes("Parent",f);
legendLabels = cell(lengthVariableParam,1);
for iter = 1:lengthVariableParam
    plot(ax,t/60,rc{iter}, "LineWidth",2)
    hold on
    legendLabels{iter} = sprintf('%s = %.0e',variableParam, variableParamArray{1}(iter));
end
set(ax,"FontSize", 28)
xlabel(ax, "Time (mins)", "FontSize", 32)
ylabel(ax, "Cavity Radius (cm)", "FontSize", 32)
xlim(ax,[t(1)/60, t(end)/60])
legend(ax, legendLabels)

% plot non-dimensional steady state solid displacement u vs non-dimesional radius
f = figure;
ax = axes("Parent",f);
legendLabels = cell(lengthVariableParam,1);
for iter = 1:lengthVariableParam
    uiter = u{iter};
    uss = uiter(:,end)/(Q/(8*c(iter)*pi));
    semilogx(ax,r/Rtot,uss, "LineWidth",2)
    hold on
    legendLabels{iter} = sprintf('%s = %.1e',variableParam, variableParamArray{1}(iter));
end
set(ax,"FontSize", 28)
xlabel(ax, "r/R_{tot}", "FontSize", 32)
ylabel(ax, "u/(Q/8c\pi)", "FontSize", 32)
xlim(ax,[1e-2, r(end)/Rtot])
legend(ax, legendLabels)

% plot steady state pressure vs radius
f = figure;
ax = axes("Parent",f);
legendLabels = cell(lengthVariableParam,1);
for iter = 1:lengthVariableParam
    piter = p{iter};
    pss = piter(:,end);
    plot(ax,r,pss/1e4, "LineWidth",2)
    hold on
    legendLabels{iter} = sprintf('%s = %.1e',variableParam, variableParamArray{1}(iter));
end
set(ax,"FontSize", 28)
xlabel(ax, "Radius (cm)", "FontSize", 32)
ylabel(ax, "Steady State Fluid Pressure (kPa)", "FontSize", 32)
xlim(ax,[0, r(end)])
legend(ax, legendLabels)

% plot pressure at r0 vs time
f = figure;
ax = axes("Parent",f);
legendLabels = cell(lengthVariableParam,1);
for iter = 1:lengthVariableParam
    piter = p{iter};
    pr0 = piter(1,:);
    plot(ax,t/60,pr0/1e4, "LineWidth",2)
    hold on
    legendLabels{iter} = sprintf('%s = %.1e',variableParam, variableParamArray{1}(iter));
end
set(ax,"FontSize", 28)
xlabel(ax, "Time (mins)", "FontSize", 32)
ylabel(ax, "Pressure at r0 (kPa)", "FontSize", 32)
xlim(ax,[0, t(end)/60])
legend(ax, legendLabels)


