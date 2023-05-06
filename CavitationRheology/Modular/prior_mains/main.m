
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

ind_t_injection_end = getIndices.getIndInjectionEnd(t, Vol, Q_needle);
ind_r0 = getIndices.getIndNeedleRadius(r, r0);

%% FINDING VELOCITIES AND VELOCITY GRADIENTS

% v_infusion = getVelocity.ethanol_velocity_v2(r, t, ind_r_cavity, ind_r0, r0, Q_needle);
% dvdr_infusion = getVelocityGradient.simple_velocity_gradient(r, ind_r0, r0, Q_needle);
% v_relaxation = getVelocity.ethanol_velocity_v2(r, t, ind_r_cavity, ind_r0, r0, 0);
% dvdr_relaxation = getVelocityGradient.simple_velocity_gradient(r, ind_r0, r0, 0);


%% COMPUTATIONS


%initializing result arrays 
tt = zeros(length(t),length(phiVec));
lambda = zeros(length(t),length(phiVec));
dlambdadt = zeros(length(t),length(phiVec));
cavity_radius = zeros(length(t), length(phiVec));
ind_r_cavity = zeros(length(t), length(phiVec));
Pc = zeros(length(t),length(phiVec));
P0 = zeros(length(t),length(phiVec));
Pinf = zeros(length(t), length(phiVec));
Pcrit = zeros(length(t), length(phiVec));
Stress = zeros(length(t), length(phiVec));
Pfield = cell(length(phiVec),1);
a = zeros(length(t), length(phiVec));
ind_a = zeros(length(t), length(phiVec));

%     Pc = getPressure.p_c(E, lambda);
%     P0 = getPressure.p_0(phi, K_oh_c, Pc, lambda, dlambdadt, Q_needle, r_needle);
%     Pinf = getPressure.p_inf(P0, mu_oh, Q_needle, L_needle, r_needle);
%     Pcrit = getPressure.p_crit(E, lambda);
%     Stress = getPressure.stress(E, lambda);

C = cell(length(phiVec),1);

%building piecewise time vectors depending on injection end
t_injection = t(1 : ind_t_injection_end);
t_relaxation = t(ind_t_injection_end + 1 : end);

% for i = 1:length(omegaVec)
%     
%     IC_oh = zeros(1,numr+1); %Initial condition
%     IC_oh(1) = y0;
%     IC_oh(2:ind_r0) = C_0;
%     
%     [t_temp, y_temp] = ode23s(@(tt,Y) getConcentration.velocity_bc(tt, Y, r, omega(i), Q_needle, r0, D_C, D_S, phi, Rtot, numr, v, dvdr, fastFind), t, IC_oh, opts1);
%     tt(:,i) = t_temp;
%     lambda(:,i) = y_temp(1,:);
% end



for i = 1:length(phiVec)
    
    %getting lambda, dlambda and cavity radius:
    %[~,lambda_temp] = ode23s(@(t,lambda) getStrainRate.velocity_bc_porosity(lambda,Q_needle,E,r0,K_oh_t,phiVec(i),theta(i),Rtot), t_injection, y0);
    %tt(:,i) = t;
    %dlambda_temp = getStrainRate.velocity_bc_porosity(lambda_temp,Q_needle,E,r0,K_oh_t,phiVec(i),Rtot);
    
    %dlambdadt(:,i) = [dlambda_temp; zeros(length(t_relaxation),1)];
    %lambda(:,i) = [lambda_temp; lambda_temp(end) * ones(length(t_relaxation),1)];
    
    %cavity_radius(:,i) = lambda(:,i) .* r0;
    %ind_r_cavity(:,i) = getIndices.getIndCavityRadius(t, r, cavity_radius(:,i));

    
    %Initial condition mass transport:
    %IC_oh = zeros(1,numr); 
    %IC_oh(1:ind_r0) = C_0;
    C_oh_inj = phiVec(i) * rho_etoh;
    IC_oh = zeros(1,numr+1); 
    IC_oh(1:ind_r0) = C_oh_inj;
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
    [t_injection, c_injection] = ode15s(@(t_injection,C_injection) getConcentration.velocity_bc_cavity_radius_concentration(t_injection, C_injection, r, D_Cvec(i), D_S, P_vec(i), theta(i), epsilon_c(i), epsilon_t, r0, Rtot, numr, Q_needle, kappa, C_oh_inj, fastFind, getVelocity, getVelocityGradient, getStrainRate, final_ind_r_cavity), t_injection, IC_oh, opts1);
    % v1.0:
    %[t_injection, c_injection] = ode15s(@(t_injection,C_injection) getConcentration.pressure_bc(t_injection, C_injection, r, D_Cvec(i), D_S, E, P_vec(i), theta(i), K_oh_t, r0, Rtot, numr, Q_needle, C_oh_inj, fastFind, getVelocity, getVelocityGradient, getPressure, final_ind_r_cavity), t_injection, IC_oh, opts1);
                                                                                    
    final_lambda = c_injection(end,end);
    final_cavity_radius = final_lambda * r0;
    final_ind_r_cavity = fastFind.binarySearchBin(r, final_cavity_radius);

    % v3.0:
    %[t_relaxation, c_relaxation] = ode15s(@(t_relaxation,C_relaxation) getConcentration.velocity_bc_cavity_radius_concentration_v2(t_relaxation, C_relaxation, r, D_Cvec(i), D_S, P_vec(i), theta(i), r0, Rtot, numr, epsilon_c(i), 0, kappa, C_oh_inj, fastFind, getVelocity, getVelocityGradient, final_ind_r_cavity), t_relaxation, c_injection(end,:), opts1);
    % v2.0:
    [t_relaxation, c_relaxation] = ode15s(@(t_relaxation,C_relaxation) getConcentration.velocity_bc_cavity_radius_concentration(t_relaxation, C_relaxation, r, D_Cvec(i), D_S, P_vec(i), theta(i), epsilon_c(i), epsilon_t, r0, Rtot, numr, 0, kappa, C_oh_inj, fastFind, getVelocity, getVelocityGradient, getStrainRate, final_ind_r_cavity), t_relaxation, c_injection(end,:), opts1);
    % v1.0:
    %[t_relaxation, c_relaxation] = ode15s(@(t_relaxation,C_relaxation) getConcentration.pressure_bc(t_relaxation, C_relaxation, r, D_Cvec(i), D_S, E, P_vec(i), theta(i), K_oh_t, r0, Rtot, numr, 0, C_oh_inj, fastFind, getVelocity, getVelocityGradient, getPressure, final_ind_r_cavity), t_relaxation, c_injection(end,:), opts1);

    lambda(:,i) = [c_injection(:,end); c_relaxation(:,end)];
    cavity_radius(:,i) = lambda(:,i) .* r0;
    tt(:,i) = [t_injection; t_relaxation];
    C{i} = [c_injection(:, 1:end-1); c_relaxation(:, 1:end-1)];
    
    toc
        
    %pressures:
    Pc(:,i) = [getPressure.p_c(E, lambda(1:length(t_injection),i)); ...
               zeros(length(t_relaxation),1)];
    P0(:,i) = [getPressure.p_0(phiVec(i), K_oh_c, Pc(1:length(t_injection),i), lambda(1:length(t_injection),i), dlambdadt(1:length(t_injection),i), Q_needle, r_needle); ...
               zeros(length(t_relaxation),1)];
    Pinf(:,i) = [getPressure.p_inf(P0(1:length(t_injection),i), mu_oh, Q_needle, L_needle, r_needle);...
                 zeros(length(t_relaxation),1)];
    Pcrit(:,i) = getPressure.p_crit(E, lambda(:,i));
    Stress(:,i) = [getPressure.stress(E, lambda(1:length(t_injection),i)); ...
                   zeros(length(t_relaxation),1)];
               
%     ind_r_cavity_out_test = zeros(length(t),1);
%     for xx = 1:length(t)
%             cavity_radius_test = ppval(cavity_radius_spline,t(xx));
%             
%             ind_r_cavity_temp = fastFind.binarySearchBin(r,cavity_radius_test);
%             if (ind_r_cavity_temp == -1)
%                 ind_r_cavity_temp = numel(r);
%             end
%             ind_r_cavity_out_test(xx) = ind_r_cavity_temp;
%     end
%         
%     figure(20+i)
%     
%     plot(t, cavity_radius(:,i),'x')
%     hold on
%     plot(t, r(ind_r_cavity_out_test),'LineWidth',2)
%     hold off
end

% vector for plotting purposes (and backing out representative velocities)
tvec=[find(t>=0*60,1) find(t>=2*60,1) find(t>=6*60,1) find(t>=30*60,1) find(t>=60*60,1) find(t>=2*60*60,1)];

% backing out ind_r_cavity, velocities, etc.
[t_val,phi_val] = size(cavity_radius);
v = cell(length(phiVec),1);

%a = (a0^3-r0^3+r(nr_cavity).^3).^(1/3)';


for j = 1:phi_val
    for i = 1:t_val
        ind_r_cavity(i,j) = fastFind.binarySearchBin(r, cavity_radius(i,j));
        a(i,j) = (a0^3-r0^3+cavity_radius(i,j).^3).^(1/3)';
        ind_a(i,j) = fastFind.binarySearchBin(r, a(i,j));
    end
    v_tvec = zeros(length(r), length(tvec));
    for i = 1:length(tvec)
        ind_r_cavity_tvec = ind_r_cavity(tvec(i),j);
        if tvec(i) > ind_t_injection_end
           v_tvec(:,i) = getVelocity.ethanol_velocity_porosity(r, ind_r_cavity_tvec, r0, theta(j), 0);
        else
           v_tvec(:,i) = getVelocity.ethanol_velocity_porosity(r, ind_r_cavity_tvec, r0, theta(j), Q_needle);
        end
    end
    v{j} = v_tvec;
end

%% plotting:


if show_cavity_radius == 1
   getPlots.plot_cavity_radius(tt, cavity_radius, phiVec, 'x_lim_end', 12*60)
end
if show_concentration_lineplots == 1
   getPlots.plot_concentration_profiles(r, tt, tvec, C, phiVec, 'group_by', 'time', 'x_lim_end', a0);
end
if show_velocity_field == 1
   time_pt = find(t>=6*60,1);
   ind_r_cavity_time_pt = ind_r_cavity(time_pt,:);
   getPlots.plot_velocity(t, time_pt, r, ind_r0, ind_r_cavity_time_pt, v, phiVec, 'x_lim_end', a0);
end
if show_mass_conservation == 1
   time_units = 'mins';
   if t(end) >= 6*60*60
       time_units = 'hrs';
   end
   getPlots.plot_mass_conservation(t, ind_a, r, ind_r0, ind_r_cavity, C, phiVec,'t_units',time_units);
end
if show_cavity_volume == 1
   getPlots.plot_cavity_volume(t, a, r, cavity_radius, phiVec); 
end
if show_ethanol_cloud_volume == 1
   time_pt = find(t>=5*60,1);
   c_threshold = 0.2 * rho_etoh;
   getPlots.plot_ethanol_cloud(t, time_pt, C, r, c_threshold, phiVec) 
end

% cavity_volume = zeros(length(t), length(phiVec));
% tumor_volume = zeros(length(t), length(phiVec));
% 
% for i = 1:phi_val
%    cavity_volume(:,i) = (4/3)*pi*cavity_radius(:,i).^3;
%    tumor_volume(:,i) = (4/3)*pi*a(:,i).^3-cavity_volume(:,i);
% end
% 
% f = figure;
% ax1 = axes('Parent', f);
% cm = colormap('copper');
% [numColors,~] = size(cm);
% cm_divisor = floor(numColors/length(phiVec));
% 
% legendLabels = cell(length(phi_val)*2,1);
% 
% for i = 1:phi_val
%     plot(ax1, tt/60, a(:,i), 'LineWidth', 4, 'color', cm(i*cm_divisor,:))
%     hold on
%     plot(ax1, tt/60, cavity_radius(:,i), 'LineWidth', 4, 'LineStyle', '--', 'color', cm(i*cm_divisor,:))
%     legendLabels{i} = sprintf('%.0f%%w/w EC, tumor',(1-phiVec(i))*100);
% end
% for label = length(phiVec)+1:2*length(phiVec)
%     legendLabels{label} = sprintf('%.0f%%w/w EC, cavity',(1-phiVec(floor(label/2)))*100);
% end
% 
% legend(ax1, legendLabels)
% set(ax1,'FontSize', 28)
% xlim([0, 6])
% ylabel('Radius (cm)', 'FontSize', 32)
% xlabel('Time (mins)', 'FontSize', 32)
% 
% f = figure;
% ax2 = axes('Parent',f);
% for i = 1:phi_val
%     plot(ax2, tt/60, cavity_volume(:,i), 'LineWidth', 4, 'color', cm(i*cm_divisor,:))
%     hold on
%     plot(ax2, tt/60, tumor_volume(:,i), 'LineWidth', 4, 'LineStyle', '--', 'color', cm(i*cm_divisor,:))
% end
% 
% legend(ax2, legendLabels)
% set(ax2,'FontSize', 28)
% xlim([0, 6])
% ylabel('Volume (cm^{3})', 'FontSize', 32)
% xlabel('Time (mins)', 'FontSize', 32)


%     figure() 
%     plot(t, ppval(cavity_radius_spline,t),'LineWidth',2)
%     hold on
%     plot(t, cavity_radius(:,i),'x')

% show_strain = 1;
% show_cavity_radius = 0;
% show_pressure_field = 1;
% show_velocity_field = 1;
% show_pc_pcrit_p0 = 0;
% show_p_inf = 0;
% show_heatmaps = 0;
% show_concentration_lineplots = 0;
% show_mass_conservation = 0;


