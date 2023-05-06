clear
clc

%% simulation params

runCavitySimulations = 0;
runConcentrationSimulations = 0;
writeResults = 0;
e_filename = 'dilatation.mat';
u_filename = 'displacement.mat';
rc_filename = 'cavity_radius.mat';
v_filename = 'velocity.mat';
dvdr_filename = 'velocity_gradient.mat';
p_filename = 'pressure.mat';
c_filename = 'concentration.mat';

plotPressureLine = 1;
plotPressure3D = 1;
plotDilatationLine = 0;
plotDisplacementLine = 0;
plotDisplacement3D = 0;
plotVelocityLine = 0;
plotVelocity3D = 0;
plotCavityRadius = 1;
plotConcentration = 1;

%% initialize params

variableParam = "K";
%K = [3.8e-11, 7.6e-11, 3.8e-10, 7.6e-10]; % m^2 / kPa = cm^2/barye
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
r = linspace(r0,Rtot,numr)';
ind_r0 = 21;
%r_inc_cavity = [linspace(0,r0-0.05,10),r];
r_mass = [linspace(0,r0-(1e-3),ind_r0-1),linspace(r0,Rtot,numr)]';
t_injection_end = Vol/Q;
t_injection = linspace(0.01,t_injection_end,600);
t_relaxation = [t_injection_end+1,t_end];

%t = 0.01:2:20000;
%t = linspace(0.01,12*60*60,400);
%rho_etoh = 0.789; % density of ethanol.
c0 = 0.5;
D_S = 5*10^(-6); %(cm^2/s)
D_C = D_S;
P = 1;

%ind_r0 = find(r_mass >= r0,1); 

c = K*(lambda+2*mu);

%% make param dictionary

paramNames = ["K","r0","Rtot","Q","lambda","mu","phi","numr","D_C","D_S","P","c0","ind_r0"];
paramIndices = 1:length(paramNames);
paramValues = {K, r0, Rtot, Q, lambda, mu, phi, numr, D_C, D_S, P, c0, ind_r0};
params = dictionary(paramNames, paramValues);
paramLabels = dictionary(paramIndices,paramNames);

%e_time = talbot_inversion_sym(@(s) getDilatationLaplace(s,r,params), t);
%u_time = talbot_inversion_sym(@(s) getDisplacementLaplace(s,r,params),

%% calculate dilatation, solid displacement, cavity radius and pressure

variableParamArray = params(variableParam);
lengthVariableParam = length(variableParamArray{1});

e = cell(lengthVariableParam,1);
u = cell(lengthVariableParam,1);
rc = cell(lengthVariableParam,1);
v = cell(lengthVariableParam,1);
dvdr = cell(lengthVariableParam,1);
p = cell(lengthVariableParam,1);
t_m = cell(lengthVariableParam,1);
c_m = cell(lengthVariableParam,1);

if runCavitySimulations == 1
    for iter = 1:lengthVariableParam
        ep = zeros(length(r),length(t_injection));
        up = zeros(length(r),length(t_injection));
        vp = zeros(length(r),length(t_injection));
        dvdrp = zeros(length(r),length(t_injection));
        pp = zeros(length(r),length(t_injection));
        paramsIterValues = getParamValues(params,paramLabels,variableParam,iter);
        paramsIter = dictionary(paramNames,paramsIterValues);

        fprintf('Continuum Mechanics: %d\n',iter)
        tic    
        for ri = 1:length(r)
            ep(ri,:) = inverseLaplace(@(s) getDilatationLaplace(s,r(ri),paramsIter), t_injection);
            up(ri,:) = inverseLaplace(@(s) getDisplacementLaplace(s,r(ri),paramsIter), t_injection);
            vp(ri,:) = inverseLaplace(@(s) getVelocityLaplace(s,r(ri),ep(ri,:),paramsIter), t_injection);
            dvdrp(ri,:) = inverseLaplace(@(s) getVelocityGradLaplace(s,r(ri),ep(ri,:),paramsIter), t_injection);
            %ep(ri,:) = matlab_ilt(@(s) getDilatationLaplace(s,r(ri),paramsIter), t, 100, 'euler');
            %up(ri,:) = matlab_ilt(@(s) getDisplacementLaplace(s,r(ri),paramsIter), t, 100, 'euler');
            %vp(ri,:) = matlab_ilt(@(s) getVelocityLaplace(s,r(ri),paramsIter), t, 100, 'euler');
            %dvdrp(ri,:) = matlab_ilt(@(s) getVelocityGradLaplace(s,r(ri),paramsIter), t, 100, 'euler');
    
            %% interpolating NaN values:
            % e
            ep(ri,1) = 0;
            ep_nan = isnan(ep(ri,:));
            er_interp = interp1(t_injection(~ep_nan),ep(ri,~ep_nan),t_injection(ep_nan),'spline','extrap');
            ep(ri,ep_nan) = er_interp;
            % u
            up(ri,1) = 0;
            up_nan = isnan(up(ri,:));
            ur_interp = interp1(t_injection(~up_nan),up(ri,~up_nan),t_injection(up_nan),'spline','extrap');
            up(ri,up_nan) = ur_interp;
            % v
            vp(ri,1) = 0;
            vp_nan = isnan(vp(ri,:));
            vr_interp = interp1(t_injection(~vp_nan),vp(ri,~vp_nan),t_injection(vp_nan),'spline','extrap');
            vp(ri,vp_nan) = vr_interp;
            % dvdr
            dvdrp(ri,1) = 0;
            dvdrp_nan = isnan(dvdrp(ri,:));
            dvdrp_interp = interp1(t_injection(~dvdrp_nan),dvdrp(ri,~dvdrp_nan),t_injection(dvdrp_nan),'spline','extrap');
            dvdrp(ri,dvdrp_nan) = dvdrp_interp;
    
            if variableParam == "mu"
                pp(ri,:) = ep(ri,:)*(lambda+2*mu(iter));
            else
                pp(ri,:) = ep(ri,:)*(lambda+2*mu);
            end
    
        end
        toc
    
        rc{iter} = r0 + up(1,:);
        e{iter} = ep;
        u{iter} = up;
        v{iter} = vp;
        p{iter} = pp;
        dvdr{iter} = dvdrp;

        %tic
        %writematrix(ep,e_filename,'Sheet',iter)
        %writematrix(up,u_filename,'Sheet',iter)
        %writematrix(rc{iter},rc_filename,'Sheet',iter)
        %writematrix(vp,v_filename,'Sheet',iter)
        %writematrix(pp,p_filename,'Sheet',iter)
        %writematrix(dvdrp,dvdr_filename,'Sheet',iter)
        %toc

    end

end

if writeResults == 1 && runCavitySimulations == 1
    tic
    fprintf('Writing Results CM: \n')
    save(e_filename,'e','t_injection')
    save(u_filename,'u','t_injection')
    save(rc_filename,'rc','t_injection')
    save(v_filename,'v','t_injection')
    save(p_filename,'p','t_injection')
    save(dvdr_filename,'dvdr','t_injection')
    toc
end

%% CONCENTRATION
% S = getSparsity(numr);
S = getSparsity(length(r_mass));
S = sparse(S);
opts1 = odeset('Vectorized','on','JPattern',S,'RelTol',1e-3,'AbsTol',1e-4,'NonNegative',true);
%opts1 = odeset('Vectorized','on','JPattern',S,'NonNegative',true);

if runConcentrationSimulations == 1
    for iter = 1:lengthVariableParam
    
        paramsIterValues = getParamValues(params,paramLabels,variableParam,iter);
        paramsIter = dictionary(paramNames,paramsIterValues);
    
        if runCavitySimulations == 0
            load(e_filename);
            load(u_filename);
            load(rc_filename);
            load(v_filename);
            load(p_filename);
            load(dvdr_filename);
            ep = e{iter};
            up = u{iter};
            rcp = rc{iter};
            vp = v{iter};
            pp = p{iter};
            dvdrp = dvdr{iter};
            %v_spline = spline(t,vp);
            %dvdr_spline = spline(t,dvdrp);
            %rc_spline = spline(t,rcp);
        else
            vp = v{iter};
            dvdrp = dvdr{iter};
            pp = p{iter};
            rcp = rc{iter};
        end
    
%         v_spline = cell(length(r),1);
%         dvdr_spline = cell(length(r),1);
%         p_spline = cell(length(r),1);
        v_spline = cell(length(r_mass),1);
        dvdr_spline = cell(length(r_mass),1);
        p_spline = cell(length(r_mass),1);
        rc_spline = spline(t_injection,rcp);

        fprintf('Spline creation: %d\n',iter)
        tic

%         for i = 1:length(r)
%             v_spline{i} = spline(t_injection,vp(i,:));
%             dvdr_spline{i} = spline(t_injection,dvdrp(i,:));
%             p_spline{i} = spline(t_injection,pp(i,:));
%         end
        for i = 1:length(r_mass)
            if i <= ind_r0
                v_spline{i} = spline(t_injection,(Q/(4*pi*r_mass(i)^2))*ones(size(vp,2),1));
                dvdr_spline{i} = spline(t_injection,(-Q/(2*pi*r_mass(i)^3))*ones(size(vp,2),1));
                p_spline{i} = spline(t_injection,pp(ind_r0,:));
            else            
                v_spline{i} = spline(t_injection,vp(i-ind_r0,:));
                dvdr_spline{i} = spline(t_injection,dvdrp(i-ind_r0,:));
                p_spline{i} = spline(t_injection,pp(i-ind_r0,:));
            end
        end

        toc
        
        fprintf('Mass Transport: %d\n',iter)
        tic
    
%         IC = zeros(1,length(r)); % initial condition
%         [t_mp,c_mp] = ode15s(@(t_p,conc_p) getConcentrationNetti(t_p, conc_p, r, v_spline, dvdr_spline, paramsIter), [t_injection(1), t_injection(end)], IC, opts1);
        IC = zeros(1,length(r_mass));
        IC(1:ind_r0) = c0;
        [t_mp,c_mp] = ode15s(@(t_p,conc_p) getConcentration_v2(t_p, conc_p, r_mass, v_spline, dvdr_spline, p_spline, paramsIter), [t_injection(1), t_injection(end)], IC, opts1);

        if t_end > t_injection(end)
            IC_relax = c_mp(end,:);
%             for i = 1:length(r)
%                 v_spline{i} = spline(t_injection,zeros(length(t_injection),1));
%                 dvdr_spline{i} = spline(t_injection,zeros(length(t_injection),1));
%                 p_spline{i} = spline(t_injection,zeros(length(t_injection),1));
%             end
            for i = 1:length(r_mass)
                v_spline{i} = spline(t_injection,zeros(length(t_injection),1));
                dvdr_spline{i} = spline(t_injection,zeros(length(t_injection),1));
                p_spline{i} = spline(t_injection,zeros(length(t_injection),1));
            end
%             [t_mp_relax,c_mp_relax] = ode15s(@(t_p,conc_p) getConcentrationNetti(t_p, conc_p, r, v_spline, dvdr_spline, paramsIter), t_relaxation, IC_relax, opts1);
            [t_mp_relax,c_mp_relax] = ode15s(@(t_p,conc_p) getConcentration_v2(t_p, conc_p, r_mass, v_spline, dvdr_spline, p_spline, paramsIter), t_relaxation, IC_relax, opts1);
            t_m{iter} = [t_mp;t_mp_relax];
            c_m{iter} = [c_mp;c_mp_relax];
        else
            t_m{iter} = t_mp;
            c_m{iter} = c_mp;
        end
    
        toc

    end
end


if writeResults == 1 && runConcentrationSimulations == 1
    fprintf('Writing Results MT: \n')
    tic
    save(c_filename,'t_m','c_m')
    toc
end

%% PLOT PRESSURE

if plotPressureLine == 1
    if runCavitySimulations == 0
        load(p_filename);
    end
    plotResults.plotPressureLine(t_injection,p,params,variableParam)
end

if plotPressure3D == 1
    if runCavitySimulations == 0
        load(p_filename);
    end
    plotResults.plotPressure3D(t_injection,r,p,params,variableParam)
end

%% PLOT DILATATION

if plotDilatationLine == 1
    if runCavitySimulations == 0
        load(e_filename);
    end
    plotResults.plotDilatationLine(t_injection,e,params,variableParam)
end

%% PLOT SOLID DISPLACEMENT

if plotDisplacementLine == 1
    if runCavitySimulations == 0
        load(u_filename);
    end
    plotResults.plotDisplacementLine(t_injection,u,params,variableParam)
end

if plotDisplacement3D == 1
    if runCavitySimulations == 0
        load(u_filename);
    end
    plotResults.plotDisplacement3D(t_injection,r,u,params,variableParam)
end

%% PLOT FLUID VELOCITY

if plotVelocityLine == 1
    if runCavitySimulations == 0
        load(v_filename);
    end
    plotResults.plotVelocityLine(t_injection,v,params,variableParam)
end

if plotVelocity3D == 1
    if runCavitySimulations == 0
        load(v_filename);
    end
    plotResults.plotVelocity3D(t_injection,r,v,params,variableParam)
end

%% PLOT CAVITY RADIUS

if plotCavityRadius == 1
    if runCavitySimulations == 0
        load(rc_filename);
    end
    plotResults.plotCavityRadius(t_injection,rc,params,variableParam)
end

%% CONCENTRATION

if plotConcentration == 1
    if runConcentrationSimulations == 0
        load(c_filename);
    end
    t_vec = [5*60, 20*60, 30*60, 60*60, 90*60];
    %t_vec = [6*60, 30*60, 60*60, 120*60];
%     plotResults.plotConcentration(t_m,t_injection,r,c_m,rc,t_vec,params,variableParam)
    plotResults.plotConcentration(t_m,t_injection,r_mass,c_m,rc,t_vec,params,variableParam)
end






