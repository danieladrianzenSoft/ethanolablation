%% TEST NUMERICAL DERIVATIVE

x = linspace(0,2*pi,100);
f = @(x,y) sin(x);
y = f(x);

dy = numericalDerivative(f, x, y, 0.1, 'x');

figure()
plot(x,y,'LineWidth',3)
hold on
plot(x,dy,'k','LineWidth',3)

%% NON-STIFF: DY/DT = 2T

tspan = [0 5];
y0 = 0;
[t,y] = ode45(@(t,y) 2*t, tspan, y0);
[tcustom, ycustom] = backwardEulerSolver(@(t,y) 2*t, y0, tspan, 0.05);

f=figure();
ax = axes('Parent',f);
plot(ax,t,y,'LineWidth',3,'color',[0.8,0.8,0.8])
hold on
plot(ax,tcustom,ycustom,'LineStyle','--','LineWidth',3,'color',[255,86,0]/255)
set(ax, 'FontSize',24)
title(ax, 'non-stiff: 2*t','FontSize',32)
legend(ax,'ODE45', 'Custom Backward Euler')

%% NON-STIFF MULTIPLE IC'S: DY/DT = -2Y + 2COS(T)SIN(2T)

yprime = @(t,y) -2*y + 2*cos(t).*sin(2*t);
y0 = -5:5; 
tspan = [0 3];
[t_mult,y_mult] = ode45(yprime,tspan,y0);
[t_mult_custom,y_mult_custom] = backwardEulerSolver(yprime,y0',tspan, 0.01);

f=figure();
ax = axes('Parent',f);
plot(ax,t_mult,y_mult,'LineWidth',3,'color',[0.8,0.8,0.8])
hold on
plot(ax,t_mult_custom,y_mult_custom,'LineWidth',3,'LineStyle','--','color',[255,86,0]/255)
set(ax, 'FontSize',24)
title(ax, 'non-stiff: y'' = -2y + 2cos(t)sin(2t) multiple ICs','FontSize',32)
legend(ax,'ODE45', '','','','','','','','','','','Custom Backward Euler')

%% NON-STIFF: VAN DER POL'S, MU = 1

[tvdp,yvdp] = ode45(@vdp1,[0 20],[2; 0]);
[tvdp_custom, ydp_custom] = backwardEulerSolver(@(t,y) vdp1(t,y), [2; 0], [0, 20], 0.005);

f=figure();
ax = axes('Parent',f);
plot(ax,tvdp,yvdp(:,1),'LineWidth',3,'color',[0.8,0.8,0.8])
hold on
plot(ax,tvdp_custom,ydp_custom(1,:),'LineWidth',3,'LineStyle','--','color',[255,86,0]/255)
set(ax, 'FontSize',24)
title(ax, 'non-stiff: van der Pol \mu = 1','FontSize',32)
legend(ax,'ODE45','Custom Backward Euler')

%% STIFF: VAN DER POL'S

dt = 0.0005;
mu = 1000;
tspan = [0, 3000];
y0 = [2,0];

%ode15s
matlabODE = tic;
[t_vdp,y_vdp] = ode15s(@(t,y) vdp(t,y,mu),tspan,y0);
matlabODEtime = toc(matlabODE);
fprintf('ODE15s: %ss \n', matlabODEtime);

%CN-explicit
%[tvdp1000_eul, ydp1000_eul] = backwardEulerSolver(@(t,y) vdp1000(t,y), [2; 0], [0, 3000], dx);
CNex = tic;
[tvdp_cn, ydp_cn] = crankNicolsonSolver(@(t,y) vdp(t,y,mu), y0', tspan, dt);
CNexTime = toc(CNex);
fprintf('Crank-Nicolson Explicit: %ss \n', CNexTime);

%CN-implicit
params = dictionary('mu',mu);
options = dictionary('tolerance',1e-3);
J = getJacobian("Vdp");
CNim = tic;
[tvdp_cni, ydp_cni] = crankNicolsonImplicitSolver(@(t,y) vdp(t,y,mu), J, y0', tspan, dt, params, helpers, options);
CNimTime = toc(CNim);
fprintf('Crank-Nicolson Implicit: %ss \n', CNimTime);

f=figure();
ax = axes('Parent',f);
plot(ax,t_vdp,y_vdp(:,1),'LineWidth',3,'color',[0.8,0.8,0.8])
hold on
%plot(ax,tvdp1000_eul,ydp1000_eul(1,:),'LineWidth',3,'LineStyle','--','color',[255,86,0]/255)
plot(ax,tvdp_cn,ydp_cn(1,:),'LineWidth',3,'LineStyle','--','color',[80,180,255]/255)
plot(ax,tvdp_cni,ydp_cni(1,:),'LineWidth',3,'LineStyle','--','color',[0.9290 0.6940 0.1250])

titleText = sprintf('stiff: van der Pol mu = %.0f', mu);
title(ax, titleText,'FontSize',32)
legend(ax,sprintf('ODE15s (%.3fs)',matlabODEtime),sprintf('Explicit C-N (%.3fs)', CNexTime), sprintf('Implicit C-N (%.3fs)', CNimTime))
set(ax, 'FontSize',24)

%% STIFF: DIFFUSION EQN

clear
clc

dx = 0.01;
dt = 1;
xSpan = [0,1];
tSpan = [0,2*60*60];
D = 1e-6;
numX = (xSpan(2)-xSpan(1)) / dx;
IC = [1*ones(numX/2,1); zeros(numX/2,1)];

%ode15s
matlabODE = tic;
[t_diff,y_diff] = ode15s(@(t,C) diffEqn1Compt(t,C,D,xSpan,dx),tSpan,IC);
matlabODETime = toc(matlabODE);
fprintf('ODE15s: %ss \n', matlabODETime);
%CN explicit
cnEx = tic; 
[t_diff_CN,y_diff_CN] = crankNicolsonSolver(@(t,C) diffEqn1Compt(t,C,D,xSpan,dx),IC,tSpan,dt);
cnExTime = toc(cnEx);
fprintf('Crank-Nicolson Explicit: %ss \n', cnExTime);
%CN implicit
solver = 'midpointImplicitSolver';
params = dictionary('D',D,'numX',numX,'dx',dx);
options = dictionary('tolerance',1e-3);
J = getJacobian("DiffEqn1Compt");
cnIm = tic;
solverName = '';
helpers = {};

JacFun = J(u_int, deltaZ, Nz, cFeed);
options = odeset( 'Jacobian' , @(t,x)JacFun(t, x, u_int, deltaZ, Nz, cFeed), 'RelTol',1e-8, 'AbsTol',1e-8);
[t,V] = ode15s(@(t,x)odeTest(t,x,u_int, deltaZ, Nz, cFeed), tspan, v0, options);

if strcmp(solver, "midpointImplicitSolver")
    [t_diff_CNI,y_diff_CNI] = midpointImplicitSolver(@(t,C) diffEqn1Compt(t,C,D,xSpan,dx),J,IC,tSpan,dt, params, helpers, options);
    solverName = "Midpoint Implicit";
elseif strcmp(solver, "crankNicolsonImplicitSolver")
    [t_diff_CNI,y_diff_CNI] = crankNicolsonImplicitSolver(@(t,C) diffEqn1Compt(t,C,D,xSpan,dx),J,IC,tSpan,dt, params, helpers, options);
    solverName = "C-N Implicit";
end
cnImTime = toc(cnIm);
fprintf('%s: %ss \n', solverName, cnImTime);

tvec = [0*60, 2*60, 30*60, 1*60*60, 1.5*60*60];


f=figure();
ax = axes('Parent',f);
x = xSpan(1):dx:xSpan(2)-dx;

for i = 1:length(tvec)
    tvec_ode15 = find(t_diff >= tvec(i),1);
    tvec_CN = find(t_diff_CN >= tvec(i),1);
    tvec_CNI = find(t_diff_CNI >= tvec(i),1);
    plot(ax,x,y_diff(tvec_ode15,:),'LineWidth',3,'color',[0.8,0.8,0.8])
    hold on
    plot(ax,x,y_diff_CN(:,tvec_CN),'LineWidth',3,'LineStyle','--','color',[80,180,255]/255)
    plot(ax,x,y_diff_CNI(:,tvec_CNI),'LineWidth',3,'LineStyle','--','color',[0.9290 0.6940 0.1250])

end

title(ax, 'stiff: diffusion eqn','FontSize',32)
legend(ax,sprintf('ODE15s (%.3fs)',matlabODETime),'','','','','','',sprintf('C-N Explicit (%.3fs)', cnExTime),'','','','','','',sprintf('%s (%.3fs)', solverName, cnImTime))
set(ax, 'FontSize',24)

%% STIFF: DIFFUSION-CONVECTION EQN IN RADIAL COORDINATES

clear
clc

R0 = 0.01;
Rtot = 2;
dr_1 = 0.0025;
dr_2 = 0.01;
rSpan = [0,Rtot];
r = [rSpan(1):dr_1:R0,R0+dr_2:dr_2:rSpan(2)];
R0_ind = find(r >= R0, 1);
numR = length(r);
dt_inj = 0.05;
dt_rel = 0.5;
Vol = 0.1; %ml
Q0 = 1/3600; %1ml/hr to ml/s.
D = 7.5*10^(-6); %(cm^2/s)

nRinit = R0_ind;
IC = zeros(numR,1);
IC(1:nRinit) = 1;

tspan_inj = [0,Vol/Q0];
tspan_rel = [Vol/Q0+dt_inj, 2*60*60];
v = [2 * Q0 / (pi * R0.^2), Q0./(4*pi*r(2:end).^2)];
%v = [0, Q0./(4*pi*r(2:end).^2)];
dvdr = [0,-Q0./(2*pi*r(2:end).^3)];

%ode15s
matlabODE = tic;
[tinj_ode15,Cinj_ode15] = ode15s(@(t,C) diffConv1Compt(t,C,D,r,v,dvdr), tspan_inj, IC);
IC_rel_ode15 = Cinj_ode15(end,:);
[trel_ode15,Crel_ode15] = ode15s(@(t,C) diffConv1Compt(t,C,D,r,zeros(1,numR),zeros(1,numR)), tspan_rel, IC_rel_ode15);
t_DIFFCON_ode15 = [tinj_ode15;trel_ode15];
C_DIFFCON_ode15 = [Cinj_ode15;Crel_ode15];

matlabODETime = toc(matlabODE);
fprintf('ODE15s: %ss \n', matlabODETime);

%implicit CN
solver = 'crankNicolsonImplicitSolver';
options = dictionary("tolerance",1e-3,"MaxIters",0);
params = dictionary('D', D, 'R0_ind', R0_ind);
helpers = {r, v, dvdr};
J = getJacobian("DiffConv1Compt");
solverName = '';
cnIm = tic;
if strcmp(solver, "midpointImplicitSolver")
    [tinj_CN,Cinj_CN] = midpointImplicitSolver(@(t,C) diffConv1Compt(t,C,D,r,v,dvdr), J, IC, tspan_inj, dt_inj, params, helpers, options);
    solverName = "Midpoint Implicit";
    IC_rel_CN = Cinj_CN(:,end);
    [trel_CN,Crel_CN] = midpointImplicitSolver(@(t,C) diffConv1Compt(t,C,D,r,zeros(1,numR),zeros(1,numR)), J, IC_rel_CN, tspan_rel, dt_rel, params, helpers, options);
elseif strcmp(solver, "crankNicolsonImplicitSolver")
    [tinj_CN,Cinj_CN] = crankNicolsonImplicitSolver(@(t,C) diffConv1Compt(t,C,D,r,v,dvdr), J, IC, tspan_inj, dt_inj, params, helpers, options);
    solverName = "C-N Implicit";
    IC_rel_CN = Cinj_CN(:,end);
    [trel_CN,Crel_CN] = crankNicolsonImplicitSolver(@(t,C) diffConv1Compt(t,C,D,r,zeros(1,numR),zeros(1,numR)), J, IC_rel_CN, tspan_rel, dt_rel, params, helpers, options);
end

cnImTime = toc(cnIm);
fprintf('%s Implicit: %ss \n', solverName, cnImTime);

t_DIFFCON_CN = [tinj_CN,trel_CN];
C_DIFFCON_CN = [Cinj_CN,Crel_CN];

tvec = [0*60, 2*60, 6*60, 30*60, 1*60*60, 1.5*60*60];
%tvec = [0*60, 1*60, 2*60, 4*60, 5*60];
f=figure();
ax = axes('Parent',f);

for i = 1:length(tvec)
    tvec_ode15 = find(t_DIFFCON_ode15 >= tvec(i),1);
    tvec_CN = find(t_DIFFCON_CN >= tvec(i),1);
    plot(ax,r,C_DIFFCON_ode15(tvec_ode15,:),'LineWidth',3,'color',[0.8,0.8,0.8])
    hold on
    plot(ax,r,C_DIFFCON_CN(:,tvec_CN),'LineWidth',3,'LineStyle','--','color',[80,180,255]/255)
end

title(ax, 'stiff: diffusion-convection in radial coords','FontSize',32)
legend(ax,sprintf('ODE15s (%.3fs)',matlabODETime),'','','','','','',sprintf('%s (%.3fs)', solverName, cnImTime))
set(ax, 'FontSize',24)


function dydt = vdp1(t,y)
%VDP1  Evaluate the van der Pol ODEs for mu = 1
%
%   See also ODE113, ODE23, ODE45.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2014 The MathWorks, Inc.

dydt = [y(2); (1-y(1)^2)*y(2)-y(1)];

end

function dydt = vdp(t,y,mu)
%VDP1000  Evaluate the van der Pol ODEs for mu = 1000.
%
%   See also ODE15S, ODE23S, ODE23T, ODE23TB.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2014 The MathWorks, Inc.

    dydt = [y(2); mu*(1-y(1)^2)*y(2)-y(1)];
end
% function dydy = dy_vdp1000(t,y)
% 
% end

function dCdt = diffEqn1Compt(t,C,D,xSpan,dx)

    numX = (xSpan(2)-xSpan(1)) / dx;
    dCdt = zeros(numX,1);
    d2Cdx2 = zeros(numX,1);

    for i = 1 %y direction BC zero flux (dcdx=0)
        dCdt(i) = 2*D*(C(i+1)-C(i))./(dx.^2);
    end
    for i = 2:numX-1 %in gel
        d2Cdx2(i) = (C(i+1)-2.*C(i)+C(i-1))./(dx.^2);
        dCdt(i) = D*d2Cdx2(i);
    end
    for i = numX %right before interface - gel
        dCdt(i) = 2*D*(C(i-1)-C(i))./(dx.^2); %BC zero flux
    end
  
end



function DCDT = diffConv1Compt(t,C,D,r,v,dvdr)

%NOTE: TO UNDERSTAND THESE EQUATIONS, LOOK UP FOR DIFFERENCE EQUATIONS: CENTRAL, BACKWARD
%AND FORWARD.

numR = length(r);
dCdt = zeros(numR,1);
dCdr = zeros(numR,1);
d2Cdr2 = zeros(numR,1);

for i = 1 %y direction BC zero flux (dcdx=0, dvdr = 0)
    dr = r(2)-r(1);
    dCdt(i) = 6*D.*(C(i+1)-C(i))./(dr.^2);
end
for i = 2:(numR-1) %in tissue
    dr = (r(i+1)-r(i));
    dCdr(i) = (C(i)-C(i-1))./dr;
    %dCdr(i) = (C(i+1)-C(i-1))./(2*dr);
    d2Cdr2(i) = (C(i+1)-2.*C(i)+C(i-1))./(dr.^2);
    dCdt(i) = ((2*D/r(i)).*dCdr(i))+(D*d2Cdr2(i))-((2/r(i)).*(v(i).*C(i)))-(C(i).*dvdr(i))-(v(i).*dCdr(i));
end
for i = numR %right before interface - stroma
    dr = r(end)-r(end-1);
    dCdt(i) = 6*D*(C(i-1)-C(i))./(dr.^2);
end

DCDT = dCdt;

end