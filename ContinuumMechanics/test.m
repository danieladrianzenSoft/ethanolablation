
%% testing

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

%% GEOMETRY

%rect = [3,4,-Rtot,Rtot,Rtot,-Rtot,-Rtot,-Rtot,Rtot,Rtot]';
outerCirc = [1,0,0,Rtot]';
circ = [1,0,0,r0]';
circ = [circ;zeros(length(outerCirc)-length(circ),1)];
gm = [outerCirc,circ];
sf = 'Cout-Cin';
ns = char('Cout','Cin');
ns = ns';
g = decsg(gm,sf,ns);

%% MODEL AND MESH

model = createpde(3); % #pdes = 3
geometryFromEdges(model,g);
figure()
pdegplot(model,"EdgeLabels","on"); 
axis equal
xlim([-1.1,1.1])

generateMesh(model,"Hmax",0.05);
figure()
pdemesh(model)
axis equal
xlim([-1.1,1.1])

%% GOVERNING EQUATIONS

% m d2u/dt2 + d du/dt - div(c*grad(u)) + au = f

%a = @(~,state) 2*hCoeff + 2*emiss*stefanBoltz*state.u.^3;
%f = 2*hCoeff*ta + 2*emiss*stefanBoltz*ta^4;
%d = thick*rho*specificHeat;
%c = thick*k;

m_coeff = zeros(3,1);
a_coeff = zeros(9,1);
a_coeff(5) = 1;
a_coeff(7) = 1;
c_coeff = zeros(6,1);
c_coeff(1) = c;
c_coeff(2) = c;
d_coeff = zeros(9,1);
d_coeff(1) = 1;
d_coeff(6) = -1;
% m_coeff = zeros(3,3);
% a_coeff = zeros(3,3);
% a_coeff(2,2) = 1;
% a_coeff(3,1) = 1;
% c_coeff = zeros(6,6);
% c_coeff(1,1) = c;
% c_coeff(2,2) = c;
% d_coeff = zeros(3,3);
% d_coeff(1,1) = 1;
% d_coeff(2,3) = -1;

%specifyCoefficients(model,"m",0,"d",d,"c",c,"a",a,"f",f);
f_func = @(location,state) f_coeff(location,state,params);
specifyCoefficients(model,'m',0,'d',d_coeff,'c',c_coeff,'a',a_coeff,'f',f_func);

%% INITIAL AND BOUNDARY CONDITIONS

setInitialConditions(model,[0;0;0]);
applyBoundaryCondition(model,"neumann","Edge",[5,6,7,8],"q",[1 0 0;1 0 0;0 0 0],"g", @(location,state) getBCg(location,state,params));
applyBoundaryCondition(model,"mixed","Edge",[1,2,3,4],"u",0,"EquationIndex",1);

%% SOLVING PDF

tlist = 0:6:30*60;
R = solvepde(model, tlist);
e = R.NodalSolution;

%% PLOTTING

timept = 0;
time_idx = find(tlist >= timept,1);
figure; 
pdeplot(model,"XYData",e(:,time_idx),"Contour","on","ColorMap","jet");
title(sprintf("Temperature In The Plate, t=%d mins",timept/60))
xlabel("X-coordinate, meters")
ylabel("Y-coordinate, meters")
axis equal

timept = 6*60;
time_idx = find(tlist >= timept,1);
figure; 
pdeplot(model,"XYData",e(:,time_idx),"Contour","on","ColorMap","jet");
title(sprintf("Temperature In The Plate, t=%d mins",timept/60))
xlabel("X-coordinate, meters")
ylabel("Y-coordinate, meters")
axis equal

timept = 30*60;
time_idx = find(tlist >= timept,1);
figure; 
pdeplot(model,"XYData",e(:,time_idx),"Contour","on","ColorMap","jet");
title(sprintf("Temperature In The Plate, t=%d mins",timept/60))
xlabel("X-coordinate, meters")
ylabel("Y-coordinate, meters")
axis equal

function g = getBCg(location,state,params)
    %mu = params("mu");
    %lam = params("lambda");
    r0 = params("r0");
    Q = params("Q");
    %BCS = [4*mu/(2*mu + lam)*(state.u(3,:)/r0);Q/(4*pi*r0^2);state.u(1,:)*r0*((2*mu+lam)/(4*mu))];
    %q = [1;1;0];
    %g = [0;Q/(4*pi*r0^2);state.u(1,:)*r0*((2*mu+lam)/(4*mu))]

    g = [0;Q/(4*pi*r0^2);state.ux(3,:)+state.uy(3,:)];
end

function f = f_coeff(location,state,params)
     f2 = (params("c")/params("phi"))*(state.ux(1,:)+state.uy(1,:));
     f3 = state.ux(3,:)+state.uy(3,:);
     %f2 = (params("c")/params("phi"))*(state.ux(1,:));
     %f3 = state.ux(3,:)+(2*state.u(3,:)/location.x);
    %f = -(v.*state.ux.*state.time + v.*state.uy.*state.time);
    f = [zeros(1,length(state.ux(1,:)));f2;f3];
end