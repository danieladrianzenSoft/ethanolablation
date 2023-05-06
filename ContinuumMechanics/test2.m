
clear
clc

% lambda = 13.16*1e4; % Netti 2003, 13.16kPa = 13.16*1e4 Barye
% mu = 6.58*1e4; % Netti 2003, 6.58kPa = 6.58*1e4 Barye
% K = 3e-10;
% c = K*(lambda+2*mu);
% phi = 0.7;
% syms e(r,t) u(r,t) v(r,t) r t
% %syms r t K lambda mu phi
% pdeeq = [diff(e,t)-c*laplacian(e); diff(u,r)+2*u/r-e; diff(u,t)-(1/(phi))*c*diff(e,r)-v];
% symCoeffs = pdeCoefficients(pdeeq,[e u v],'Symbolic',true)

% m d2u/dt2 + d du/dt - div(c*grad(u)) + au = f

%% PARAMETERS

% k = 400; % thermal conductivity of copper, W/(m-K)
% rho = 8960; % density of copper, kg/m^3
% specificHeat = 386; % specific heat of copper, J/(kg-K)
% thick = .01; % plate thickness in meters
% stefanBoltz = 5.670373e-8; % Stefan-Boltzmann constant, W/(m^2-K^4)
% hCoeff = 1; % Convection coefficient, W/(m^2-K)
% % The ambient temperature is assumed to be 300 degrees-Kelvin.
% ta = 300;
% emiss = .5; % emissivity of the plate surface
Q = 1/3600;
Vol = 0.50; %total injection volume -> CHANGE TO 0.5ML - 2.5ML for humans
K = 3e-10;
lambda = 13.16*1e4; % Netti 2003, 13.16kPa = 13.16*1e4 Barye
mu = 6.58*1e4; % Netti 2003, 6.58kPa = 6.58*1e4 Barye
c = K*(lambda+2*mu);
phi = 0.7;
r0 = 0.03;
Rtot = 1;
c0 = 0.5;
D_S = 5*10^(-6); %(cm^2/s)
D_C = D_S;
P = 1;

paramNames = ["K","r0","Rtot","Q","lambda","mu","c","phi","D_C","D_S","P","c0"];
paramValues = [K, r0, Rtot, Q, lambda, mu, c, phi, D_C, D_S, P, c0];
params = dictionary(paramNames, paramValues);
%paramIndices = 1:length(paramNames);
%paramLabels = dictionary(paramIndices,paramNames);


function dedt = getDilatation(t,e,r,params)

    dr = r(2)-r(1);
    dedr = zeros(length(r),1);
    d2edr2 = zeros(length(r),1);
    dedt = zeros(length(r),1);

    c = params('c');

    % r = r0
    dedr(1,:) = (e(2,:)-e(1,:))/dr;
    d2edr2(1,:) = e(2,:)-2*e()
    dedt(1,:) = c*d2edr2(1,:) + (2./r(1,:)).*dedr(1,:);
    % r
    dedt()

    if i == 1
        
    end

end

function [t,y] = crankNicolsonImplicitSolver(f, y0, tSpan, dt, options)

    if isKey(options,"tolerance") == 0
       options("tolerance") = 1e-6;
    end
    if isKey(options,"maxIters") == 0
       options("maxIters") = 200; 
    end
    
    numT = floor((tSpan(2)-tSpan(1)) / dt) + 1;

    t = [tSpan(1); zeros(numT,1)];
    y = [y0, zeros(length(y0), numT)];
    
    for i = 1:numT
        t(i+1) = t(i) + dt;
        y_guess = y(:,i) + dt * f(t(i), y(:,i));
        y(:,i+1) = newtonRaphson(f, t(i), t(i+1), y(:,i), y_guess, dt, options("tolerance"), options("maxIters"));
    end
end

function Y = newtonRaphson(f, tcurr, tnext, ycurr, Y, dt, tolerance, maxIters)
    error = 1;
    iter = 1;
    while (error >= tolerance && iter <= maxIters)
        [F, jcurr, jnext] = crankNicolsonStep(tcurr, tnext, ycurr, Y, dt, f);
        diff = ((dt/2)*(jcurr + jnext)-eye(numel(Y)))\F;
        Y = Y - diff;
        error = norm(diff, inf) / norm(Y,inf);
        iter = iter + 1;
    end
end

function [F, jcurr, jnext] = crankNicolsonStep(tcurr, tnext, ycurr, ynext, dt, f)
    [fnext, jnext] = f(tnext, ynext);
    [fcurr, jcurr] = f(tcurr, ycurr); 
    F = ycurr + (dt / 2) * (fnext + fcurr) - ynext;
end