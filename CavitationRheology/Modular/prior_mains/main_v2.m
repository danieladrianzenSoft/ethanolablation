% uses cell arrays to account for different lengths of the time vector
% which results from different values of dt selected by the adaptive
% approach of ode15s. 


%% PARAMETERS: initializeParams.m

initializeParams

%% IMPORTING REQUIRED FUNCTION LIBRARIES

fastFind = fastFindLib;
getIndices = getIndicesLib;
getSparsity = getSparsityLib;
getStrainRate = getStrainRateLib;
getPressure = getPressureLib;
getPressureGradient = getPressureGradientLib;
getVelocity = getVelocityLib;
getVelocityGradient = getVelocityGradientLib;
getConcentration = getConcentrationLib;
getPlots = getPlotsLib;

%% indices:

%ind_t_injection_end = getIndices.getIndInjectionEnd(t, Vol, Q_needle);
ind_r_needle = getIndices.getIndNeedleRadius(r, r_needle);

%% FINDING VELOCITIES AND VELOCITY GRADIENTS

% v_infusion = getVelocity.ethanol_velocity_v2(r, t, ind_r_cavity, ind_r0, r0, Q_needle);
% dvdr_infusion = getVelocityGradient.simple_velocity_gradient(r, ind_r0, r0, Q_needle);
% v_relaxation = getVelocity.ethanol_velocity_v2(r, t, ind_r_cavity, ind_r0, r0, 0);
% dvdr_relaxation = getVelocityGradient.simple_velocity_gradient(r, ind_r0, r0, 0);


%% COMPUTATIONS


%initializing result arrays 
tt = cell(length(phiVec),1);
lambda = cell(length(phiVec),1);
dlambdadt = cell(length(phiVec),1);
cavity_radius = cell(length(phiVec),1);
ind_r_cavity = cell(length(phiVec),1);
Pc = cell(length(phiVec),1);
P0 = cell(length(phiVec),1);
Pinf = cell(length(phiVec),1);
Pcrit = cell(length(phiVec),1);
Stress = cell(length(phiVec),1);
Pfield = cell(length(phiVec),1);
a = cell(length(phiVec),1);
ind_a = cell(length(phiVec),1);


C = cell(length(phiVec),1);

%building piecewise time vectors depending on injection end
%t_injection = t(1 : ind_t_injection_end);
%t_relaxation = t(ind_t_injection_end + 1 : end);

%building piecewise time spans depending on injection end
t_injection_span = [t_injection_start, t_injection_end];
t_relaxation_span = [t_injection_end + dt_injection, t_end];
tt_injection_relaxation = cell(length(phiVec),2);


for i = 1:length(phiVec)

    %Initial condition mass transport:
    C_oh_inj = phiVec(i) * rho_etoh;
    IC_oh = zeros(1,numr+1); 
    IC_oh(1:ind_r_needle) = C_oh_inj;
    IC_oh(end) = 1;
        
    %fitting spline to use in mass transport:
    %cavity_radius_spline = spline(t, cavity_radius(:,i));

    %mass transport: 
    tic
    % infusion
    final_ind_r_cavity = 1;
    % v3.0:
    %[t_injection, c_injection] = ode15s(@(t_injection,C_injection) getConcentration.velocity_bc_cavity_radius_concentration_v2(t_injection, C_injection, r, D_Cvec(i), D_S, P_vec(i), theta(i), r0, Rtot, numr, epsilon_c(i), Q_needle, kappa, C_oh_inj, fastFind, getVelocity, getVelocityGradient, final_ind_r_cavity), t_injection, IC_oh, opts1);
    % v2.0:
    [t_injection, c_injection] = ode15s(@(t_injection,C_injection) getConcentration.velocity_bc_cavity_radius_concentration(t_injection, C_injection, r, D_C_vec(i), D_S, P_vec(i), theta_vec(i), epsilon_c_vec(i), epsilon_t, r_needle, Rtot, numr, Q_needle, kappa, C_oh_inj, fastFind, getVelocity, getVelocityGradient, getStrainRate, final_ind_r_cavity), t_injection_span, IC_oh, opts1);
    % v1.0:
    %[t_injection, c_injection] = ode15s(@(t_injection,C_injection) getConcentration.pressure_bc(t_injection, C_injection, r, D_Cvec(i), D_S, E, P_vec(i), theta(i), K_oh_t, r0, Rtot, numr, Q_needle, C_oh_inj, fastFind, getVelocity, getVelocityGradient, getPressure, final_ind_r_cavity), t_injection, IC_oh, opts1);
                                                                                    
    final_lambda = c_injection(end,end);
    final_cavity_radius = final_lambda * r_needle;
    final_ind_r_cavity = fastFind.binarySearchBin(r, final_cavity_radius);
    
    IC_relaxation = c_injection(end,:);

    % v3.0:
    %[t_relaxation, c_relaxation] = ode15s(@(t_relaxation,C_relaxation) getConcentration.velocity_bc_cavity_radius_concentration_v2(t_relaxation, C_relaxation, r, D_Cvec(i), D_S, P_vec(i), theta(i), r0, Rtot, numr, epsilon_c(i), 0, kappa, C_oh_inj, fastFind, getVelocity, getVelocityGradient, final_ind_r_cavity), t_relaxation, c_injection(end,:), opts1);
    % v2.0:
    [t_relaxation, c_relaxation] = ode15s(@(t_relaxation,C_relaxation) getConcentration.velocity_bc_cavity_radius_concentration(t_relaxation, C_relaxation, r, D_C_vec(i), D_S, P_vec(i), theta_vec(i), epsilon_c_vec(i), epsilon_t, r_needle, Rtot, numr, 0, kappa, C_oh_inj, fastFind, getVelocity, getVelocityGradient, getStrainRate, final_ind_r_cavity), t_relaxation_span, IC_relaxation, opts1);
    % v1.0:
    %[t_relaxation, c_relaxation] = ode15s(@(t_relaxation,C_relaxation) getConcentration.pressure_bc(t_relaxation, C_relaxation, r, D_Cvec(i), D_S, E, P_vec(i), theta(i), K_oh_t, r0, Rtot, numr, 0, C_oh_inj, fastFind, getVelocity, getVelocityGradient, getPressure, final_ind_r_cavity), t_relaxation, c_injection(end,:), opts1);

    lambda{i} = [c_injection(:,end); c_relaxation(:,end)];
    cavity_radius{i} = lambda{i} .* r_needle;
    tt_injection_relaxation{i,1} = t_injection;
    tt_injection_relaxation{i,2} = t_relaxation;
    tt{i} = [t_injection; t_relaxation];
    C{i} = [c_injection(:, 1:end-1); c_relaxation(:, 1:end-1)];
    
    toc
        
    %pressures:

               
end

% PRESSURES AND STRESS
for i = 1:length(phiVec)
    C_vec = C{i};
    lambda_vec = lambda{i};
    tt_vec = tt{i};
    cavity_radius_vec = cavity_radius{i};
    dlambdadt_vec = zeros(length(tt_vec),1);
    tt_injection_vec = tt_injection_relaxation{i,1};
    tt_relaxation_vec = tt_injection_relaxation{i,2};
    
    Pc_vec = [getPressure.p_c(E, lambda_vec(1:length(tt_injection_vec))); ...
              zeros(length(tt_relaxation_vec),1)];
    
    for j = 1:length(tt_vec)
        ind_r_cavity_value = fastFind.binarySearchBin(r, cavity_radius_vec(j));
        dlambdadt_vec(j) = getStrainRate.concentration_strain_rate_porosities(r, C_vec(j,:)', C_vec(j,ind_r_cavity_value)', ind_r_cavity_value, lambda_vec(j), r_needle, kappa, Q_needle, C_oh_inj, epsilon_c_vec(i), epsilon_t);
    end

    Pc{i} = Pc_vec;
    P0_vec = [getPressure.p_0(phiVec(i), K_oh_c, Pc_vec(1:length(tt_injection_vec)), lambda_vec(1:length(tt_injection_vec)), dlambdadt_vec(1:length(tt_injection_vec)), Q_needle, r_needle); ...
               zeros(length(tt_relaxation_vec),1)];
    P0{i} = P0_vec;
    Pinf{i} = [getPressure.p_inf(P0_vec(1:length(tt_injection_vec)), mu_oh, Q_needle, L_needle, r_needle);...
                 zeros(length(tt_relaxation_vec),1)];
    Pcrit{i} = getPressure.p_crit(E, lambda_vec);
    Stress{i} = [getPressure.stress(E, lambda_vec(1:length(tt_injection_vec))); ...
                   zeros(length(tt_relaxation_vec),1)];
    
end

% vector for plotting purposes (and backing out representative velocities)
%tvec=[find(t>=0*60,1) find(t>=2*60,1) find(t>=6*60,1) find(t>=30*60,1) find(t>=60*60,1) find(t>=2*60*60,1)];
tvec = [0*60, 2*60, 6*60, 30*60, 60*60, 2*60*60];
tvec_indices = cell(length(phiVec), 1);

% backing out ind_r_cavity, velocities, %tvec_indices etc.
v = cell(length(phiVec),1);


for j = 1:length(phiVec)
    t_nums = numel(tt{j});
    tt_vec = tt{j};
    cavity_radius_phi = cavity_radius{j};
    ind_r_cavity_phi = zeros(t_nums, 1);
    a_phi = zeros(t_nums, 1);
    ind_a_phi = zeros(t_nums, 1);
    for i = 1:t_nums
        ind_r_cavity_phi(i) = fastFind.binarySearchBin(r, cavity_radius_phi(i));
        a_phi(i) = (a0^3-r_needle^3+cavity_radius_phi(i).^3).^(1/3)';
        ind_a_phi(i) = fastFind.binarySearchBin(r, a_phi(i));
    end
    %cavity_radius{j} = cavity_radius_phi;
    ind_r_cavity{j} = ind_r_cavity_phi;
    a{j} = a_phi;
    ind_a{j} = ind_a_phi;
    
    v_tvec = zeros(length(r), length(tvec));
    tvec_indices_phi = zeros(length(tvec), 1);

    for i = 1:length(tvec)
        tvec_indices_phi(i) = find(tt{j} >= tvec(i),1);
        ind_r_cavity_tvec = ind_r_cavity_phi(tvec_indices_phi(i));
        if tvec(i) > t_injection_end
           v_tvec(:,i) = getVelocity.ethanol_velocity_porosity(r, ind_r_cavity_tvec, r_needle, theta_vec(j), 0);
        else
           v_tvec(:,i) = getVelocity.ethanol_velocity_porosity(r, ind_r_cavity_tvec, r_needle, theta_vec(j), Q_needle);
        end
    end
    tvec_indices{j} = tvec_indices_phi;
    v{j} = v_tvec;
end

%% plotting:


if show_cavity_radius == 1
   getPlots.plot_cavity_radius_cell(tt, cavity_radius, phiVec, 'x_lim_end', 12*60)
end
if show_concentration_lineplots == 1
   getPlots.plot_concentration_profiles_cell(r, tt, tvec, tvec_indices, C, phiVec, 'group_by', 'time', 'x_lim_end', a0);
end
% if show_velocity_field == 1
%    time_pts = 6*60;
%    tvec_time_pts_inds = cell(length(phiVec),1);
%    ind_r_cavity_time_phi = cell(length(phiVec),1);
%    for i = 1:length(phiVec)
%        ind_r_cavity_phi = ind_r_cavity{i};
%        ind_r_cavity_time_vec = zeros(length(time_pts),1);
%        tvec_time_pts_inds_phi = zeros(length(time_pts),1);
%        time_phi = tt{i};
%        for j = 1:length(time_pts)
%            time_phi_idx = find(time_phi >= time_pts(j),1);
%            tvec_time_pts_inds_phi(j) = time_phi_idx;
%            ind_r_cavity_time_vec(j) = ind_r_cavity_phi(time_phi_idx);
%        end
%        tvec_time_pts_inds{i} = tvec_time_pts_inds_phi;
%        ind_r_cavity_time_phi{i} = ind_r_cavity_time_vec;
%    end
%    %ind_r_cavity_time_pt = ind_r_cavity(time_pt,:);
%    %getPlots.plot_velocity(t, time_pts, r, ind_r0, ind_r_cavity_time_pt, v, phiVec, 'x_lim_end', a0);
%    getPlots.plot_velocity_cell(t, time_pts, tvec_time_pts_inds, r, ind_r0, ind_r_cavity_time_phi, v, phiVec, 'x_lim_end', a0);
% 
% end
if show_mass_conservation == 1
   time_units = 'mins';
   if t_relaxation(end) >= 6*60*60
       time_units = 'hrs';
   end
   getPlots.plot_mass_conservation_cell(tt, ind_a, r, ind_r_needle, ind_r_cavity, C, phiVec,'t_units',time_units);
end
if show_cavity_volume == 1
   getPlots.plot_cavity_volume_cell(tt, a, r, cavity_radius, phiVec); 
end
if show_ethanol_cloud_volume == 1
   time_pts = 5*60;
   t_vec_time_pts = cell(length(phiVec),1);
   for i = 1:length(phiVec)
       t_vec_time_phi = zeros(length(time_pts),1);
       t_phi = tt{i};
       for j = 1:length(time_pts)
           t_vec_time_phi(j) = find(t_phi >= time_pts(j),1);
       end
       t_vec_time_pts{i} = t_vec_time_phi;

   end
   c_threshold = 0.2 * rho_etoh;
   getPlots.plot_ethanol_cloud_cell(tt, time_pts, t_vec_time_pts, C, r, c_threshold, phiVec) 
end
if show_pc_pcrit == 1
   time_units = 'mins';
   getPlots.plot_pc(tt, Pc, Pcrit, phiVec, Vol, Q_needle, 't_units', time_units) 
end


