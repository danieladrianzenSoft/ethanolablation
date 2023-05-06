classdef getStrainRateLib
    methods(Static)
        function dLdt = strain_rate_newtonian_fluid_plus_solid(r, C, C_intf, ind_r_cavity, lambda, r0, kappa, Q_needle, C_oh_inj, epsilon_c, epsilon_t)
            % continuous velocity B.C, using concentration in mass
            % conservation and porosities.
            sum_mass_cavity = 4 * pi * trapz(r(1:ind_r_cavity), (r(1:ind_r_cavity)'.^2) .* C(1:ind_r_cavity,:), 1);
            vol_cavity = (4 / 3) * pi * r(ind_r_cavity).^3;
            C_avg_cavity = sum_mass_cavity ./ vol_cavity;
            dLdt = (Q_needle * (C_oh_inj - 4 * epsilon_t .* C_intf)) ./ (4 * pi * epsilon_c * r0^3 * (kappa + 1) * lambda.^2 .* C_avg_cavity); 
        end
        function dLdt = strain_rate_newtonian_fluid(r, C, C_intf, ind_r_cavity, lambda, r0, kappa, mu, Q_needle, C_oh_inj, theta, epsilon_c, epsilon_t, Pc)
            % continuous velocity B.C, using concentration in mass
            % conservation and porosities.
            sum_mass_cavity = 4 * pi * trapz(r(1:ind_r_cavity), (r(1:ind_r_cavity)'.^2) .* C(1:ind_r_cavity,:), 1);
            vol_cavity = (4 / 3) * pi * r(ind_r_cavity).^3;
            C_avg_cavity = sum_mass_cavity ./ vol_cavity;
            CC = ((3 * mu * lambda.^2 * r0) ./ (3*lambda.^2)) * ((Q_needle / (pi * r0^2)) - ((Pc * r0) / (3 * mu)));
            dLdt = (C_oh_inj * Q_needle) ./ (4 * pi * r0^3 * epsilon_c * lambda.^2 * C_avg_cavity) - (C_intf * epsilon_t * (lambda.^2 * r0^2 * Pc + 2 * CC))./(3 * mu * epsilon_c * r0^2 * lambda * C_avg_cavity); 
        end
        function dLdt = concentration_strain_rate_source(r, C, ind_r_cavity, lambda, r0, kappa, Q_needle, C_oh_inj, theta, epsilon_c)
            % v3.0
            sum_mass_cavity = 4 * pi * trapz(r(1:ind_Y_r_cavity), (r(1:ind_r_cavity)'.^2) .* C(1:ind_r_cavity,:), 1);
            vol_cavity = (4 / 3) * pi * r(ind_r_cavity).^3;
            C_avg_cavity = sum_mass_cavity ./ vol_cavity;
            dLdt = Q_needle * (C_oh_inj - (theta ./ epsilon_c) .* (lambda.^3 + epsilon_c - 1) .* C(ind_r_cavity+1,:)) ./ (4 * pi * r0^3 * (kappa + 1) * lambda.^2 .* C_avg_cavity);
        end
        function [dLdt, J] = concentration_strain_rate_porosities_jacobian(r, C, C_intf, ind_r_cavity, lambda, r0, kappa, Q_needle, C_oh_inj, epsilon_c, epsilon_t)
            % v2.1
            sum_mass_cavity = 4 * pi * trapz(r(1:ind_r_cavity), (r(1:ind_r_cavity)'.^2) .* C(1:ind_r_cavity,:), 1);
            vol_cavity = (4 / 3) * pi * r(ind_r_cavity).^3;
            C_avg_cavity = sum_mass_cavity ./ vol_cavity;
            %dLdt = Q_needle * (C_oh_inj - epsilon_t .* C(ind_r_cavity+1,:)) ./ (4 * pi * r0^3 * (kappa + 1) * lambda.^2 .* C_avg_cavity * epsilon_c);
            dLdt = Q_needle * (C_oh_inj - epsilon_t .* C_intf) ./ (4 * pi * r0^3 * (kappa + 1) * lambda.^2 .* C_avg_cavity * epsilon_c);
            J = -Q_needle * (C_oh_inj - epsilon_t .* C_intf) ./ (2 * pi * r0^3 * (kappa + 1) * lambda.^3 .* C_avg_cavity * epsilon_c);
        end
        function dLdt = concentration_strain_rate_porosities(r, C, C_intf, ind_r_cavity, lambda, r0, kappa, Q_needle, C_oh_inj, epsilon_c, epsilon_t)
            % v2.1
            sum_mass_cavity = 4 * pi * trapz(r(1:ind_r_cavity), (r(1:ind_r_cavity)'.^2) .* C(1:ind_r_cavity,:), 1);
            vol_cavity = (4 / 3) * pi * r(ind_r_cavity).^3;
            C_avg_cavity = sum_mass_cavity ./ vol_cavity;
            %dLdt = Q_needle * (C_oh_inj - epsilon_t .* C(ind_r_cavity+1,:)) ./ (4 * pi * r0^3 * (kappa + 1) * lambda.^2 .* C_avg_cavity * epsilon_c);
            dLdt = Q_needle * (C_oh_inj - epsilon_t .* C_intf) ./ (4 * pi * r0^3 * (kappa + 1) * lambda.^2 .* C_avg_cavity * epsilon_c);
        end
        function dLdt = concentration_strain_rate(r, C, ind_r_cavity, lambda, r0, kappa, Q_needle, C_oh_inj, theta)
            % v2.0
            sum_mass_cavity = 4 * pi * trapz(r(1:ind_r_cavity), (r(1:ind_r_cavity)'.^2) .* C(1:ind_r_cavity,:), 1);
            vol_cavity = (4 / 3) * pi * r(ind_r_cavity).^2;
            C_avg_cavity = sum_mass_cavity ./ vol_cavity;
            dLdt = Q_needle * (C_oh_inj - theta .* C(ind_r_cavity+1,:)) ./ (4 * pi * r0^3 * (kappa + 1) * lambda.^2 .* C_avg_cavity);
        end
        function dLdt = velocity_bc_porosity(lambda,Q_needle,E,r0,k_T,phi,theta,a)
            % CONTINUOUS VELOCITY B.C. AT CAVITY BOUNDARY 
            dLdt = ((Q_needle*(1-phi*theta))./(4*pi*(r0^3)*(lambda.^2)));
        end
        function dLdt = velocity_bc(lambda,Q_needle,E,r0,k_T,phi,a)
            % CONTINUOUS VELOCITY B.C. AT CAVITY BOUNDARY 
            dLdt = ((Q_needle*(1-phi))./(4*pi*(r0^3)*(lambda.^2)));
        end
        function dLdt = pressure_bc(lambda,Q_needle,r0,k_T,phi,a,p_c)
            % CONTINUOUS PRESSURE B.C. AT CAVITY BOUNDARY 
            dLdt = (Q_needle./(4*pi*(r0^3)*(lambda.^2)))-((k_T*p_c*phi./((r0^3).*(lambda.^2))).*(a./((a./(lambda*r0))-1)));
        end

    end
end


