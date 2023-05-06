clear
clc

%% initialize params

%variableParam = "K";
%K = [3.8e-11, 7.6e-11, 3.8e-10, 7.6e-10]; % m^2 / (kPa s) = cm^2 /(barye s)
%K = 1.5e-11;
K = 3e-10;
r0 = 0.03; % Netti 2003, 0.035cm
%r0 = round(0.311/(2*10),3); %nominal inner diameter of 24 gauge needle is 0.2032mm. Here we change to cm. 
Rtot = 1; % Netti 2003, cm
%Q = (0.1/1e3)/60; % 0.1uL/min = 0.1mm^3/min Netti 2003
Q = 1/3600;
Vol = 0.50; %total injection volume -> CHANGE TO 0.5ML - 2.5ML for humans
lambda = 13.16*1e4; % Netti 2003, 13.16kPa = 13.16*1e4 Barye
%mu = [6.58*1e4,1*1e5]; % Netti 2003, 6.58kPa = 6.58*1e4 Barye
mu = 6.58*1e4;
phi = 0.2; % Netti 2003, Porosity, dimensionless
numr = 400;
t_end = 2*60*60;
t_injection_end = Vol/Q;
t_injection = linspace(0.01,t_injection_end,600);
t_injection_span = [0,t_injection(end)];
t_relaxation_span = [t_injection_end+1,t_end];

%t = 0.01:2:20000;
%t = linspace(0.01,12*60*60,400);
%rho_etoh = 0.789; % density of ethanol.
c0 = 0.5;
D_S = 5*10^(-6); %(cm^2/s)
D_C = D_S;
P = 1;
%gamma = (0.01/((2*pi)^2));
%gamma = 1;
gamma = 0.01;

rc = 0.3;
rs = 0.2;

% values
Nx=100;
Ny=100;
x = linspace(-1,1,Nx);
y = linspace(-1,1,Ny);
dx = x(2)-x(1);
dy = y(2)-y(1);
dr = sqrt(dx^2+dy^2);
u0=zeros(Nx,Ny);
i0=0;
j0=0;

% %% DIFFUSION MASS TRANSPORT
% 
% numr = Nx;
% r = linspace(0,1,numr);
% 
% %rc_ind = find(r>=rc,1);
% %IC = zeros(length(r),1);
% %IC(1:rc_ind) = c0;
% 
% S = getSparsity(length(r));
% S = sparse(S);
% opts1 = odeset('JPattern',S,'Vectorized','on','RelTol',1e-3,'AbsTol',1e-4);
% 
% paramNames = ["K","r0","Rtot","Q","lambda","mu","phi","numr","D_C","D_S","P","c0"];
% paramIndices = 1:length(paramNames);
% paramValues = {K, r0, Rtot, Q, lambda, mu, phi, numr, D_C, D_S, P, c0};
% params = dictionary(paramNames, paramValues);
% paramLabels = dictionary(paramIndices,paramNames);
% variableParam = "K";
% 
% v_spline = cell(length(r),1);
% dvdr_spline = cell(length(r),1);
% 
% %IC_relax = zeros(length(r),1);
% %IC_relax(1:rc_ind) = c0;
% 
% load('concentration.mat')
% c_mp = c_m{1};
% t_mp = t_m{1};
% t_ind = find(t_mp>=t_injection_end,1);
% buffer = c0*ones(size(c_mp,1),10);
% c_mp = [buffer,c_mp];
% r_mp = [linspace(0,r0-0.05,10),linspace(r0,Rtot,400)];
% c_mp_spline = spline(r_mp,c_mp(t_ind,:));
% IC_relax = ppval(r,c_mp_spline);
% ind_r0 = find(r >= r0,1);
% 
% IC_relax(1:ind_r0) = c0;
% 
% %IC_relax = c_mp(end,:);
% 
% for i = 1:length(r)
%     v_spline{i} = spline(t_injection,zeros(length(t_injection),1));
%     dvdr_spline{i} = spline(t_injection,zeros(length(t_injection),1));
% end
% paramsIterValues = getParamValues(params,paramLabels,variableParam,1);
% paramsIter = dictionary(paramNames,paramsIterValues);
% 
% [t_mp_relax,c_mp_relax] = ode15s(@(t_p,conc_p) getConcentrationDiffusion(t_p, conc_p, r, v_spline, dvdr_spline, paramsIter), t_relaxation, IC_relax, opts1);
% %t_m{iter} = [t_mp;t_mp_relax];
% %c_m{iter} = [c_mp;c_mp_relax];
% 
% % tic
% % [t_relax,ut] = ode15s(@(t,u) getConcentrationDiffusion(t, u, r, D_S, c0, rc_ind, numr), t_relaxation, IC, opts1);
% % toc
% 
% 
% t_vec = [5*60, 20*60, 30*60, 60*60, 90*60];
% %t_vec = [6*60, 30*60, 60*60, 120*60];
% 
% plotResults.plotConcentration(t_mp_relax,r,c_mp_relax,t_vec,params,variableParam)


%% CAHN HILLIARD

%load('cavity_radius.mat');
%cavity_radius = rc{1};

temp_t = [0,2*60*60];
temp_rc = 0.3*ones(2,1);
rc_spline = spline(temp_t,temp_rc);

%x = 1:2:10;
%y = 1:5;

%[xx,yy] = meshgrid(x,y);
%domain_x = x(2:length(x)-1);
%domain_y = y(2:length(y)-1);
%domain_x_rep = repmat(domain_x,1,length(domain_y));
%domain_y_rep = repelem(domain_y,length(domain_x));

%domain_x_mesh = reshape(xx(2:end-1,2:end-1)',[],1)';
%domain_y_mesh = reshape(yy(2:end-1,2:end-1)',[],1)';
% checking if rc spline worked
%temp_t_2 = 0:1:2*60;
%temp_rc_2 = ppval(rc_spline,temp_t_2);
%figure()
%plot(temp_t_2,temp_rc_2);

%%

% making grid
[xx,yy]=meshgrid(x,y);
% getting mask of cavity (including rc boundary)
c_mask = (xx-i0).^2+(yy-j0).^2<=rc^2;

% getting mask of rc boundary
u0(c_mask)=1;
rc_perimeter = bwperim(u0);

% getting mask of radii inside and outside rc_perimeter, and then obtaining
% column indices
[oPerim,iPerim] = getAdjacentPerimeters(u0,rc_perimeter);
oPerim_inds = mask2ind(logical(oPerim));
iPerim_inds = mask2ind(logical(iPerim));
perim = mask2ind(rc_perimeter);

% setting inditial condition: 0 everywhere outside cavity, 1 on radius and
% also on inside perimeter
u0(~rc_perimeter)=0;
u0(rc_perimeter)=1;
u0(iPerim_inds)=1;

% getting indices of simulation domain (inside cavity, excluding
% rc_perimeter).
[c_mask_x,c_mask_y] = find(c_mask==0);
domain_inds = sub2ind([length(y),length(x)],c_mask_x,c_mask_y);

% running cahn-hilliard simulation

%D_CH = 1e-8;
D_CH = 1e-6;

%spinodal_decomp(D_CH,gamma,x,y,u0,rc_perimeter,logical(oPerim),logical(iPerim),c_inds);
%spinodal_decomp(D_CH,gamma,x,y,u0,perim_inds,oPerim_inds,iPerim_inds,c_inds);
spinodal_decomp(D_CH,gamma,x,y,t_injection_span,rc_spline,u0,rc_perimeter,c_mask,"RunNewSimulation","true");
%spinodal_decomp(D_CH,gamma);




