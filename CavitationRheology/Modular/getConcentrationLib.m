classdef getConcentrationLib
    methods(Static)
        function [F, J] = velocity_bc_cavity_radius_concentration_with_jacobian(t, Y, r, params, helpers)
            % HERE WE SOLVE FOR THE STRAIN LAMBDA AND THE CONCENTRATION
            % SIMULTANEOUSLY
            dYdt = zeros(numel(r)+1,size(Y,2));
            dYdr = zeros(numel(r)+1,size(Y,2));
            d2Ydr2 = zeros(numel(r)+1,size(Y,2));
            
            j = zeros(numel(r)+1, numel(r)+1);

            final_ind_r_cavity = params("final_ind_r_cavity");
            r_needle = params("r_needle");
            Q_needle = params("Q_needle");
            numr = params("numr");
            D_C = params("D_C");
            D_S = params("D_S");
            P = params("P");
            Rtot = params("Rtot");
            kappa = params("kappa");
            epsilon_t = params("epsilon_t");
            epsilon_c = params("epsilon_c");
            C_oh_inj = params("C_oh_inj");

            if final_ind_r_cavity == 1

                ind_r_cavity = helpers{1}(r, Y(end,end) * r_needle);

                if (ind_r_cavity == -1)
                    ind_r_cavity = final_ind_r_cavity;
                end
       
            else
                ind_r_cavity = final_ind_r_cavity;
                
            end
            
            ind_Y_r_cavity = ind_r_cavity;

            v = helpers{2}(r, ind_r_cavity, r_needle, 1, Q_needle);
            dvdr = helpers{3}(r, ind_r_cavity, 1, Q_needle);
            
            %Concentration at cavity/stroma interface
            C_intf_CSa = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_S./P)+D_C); %Right before interface
            C_intf_CSb = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_C.*P)+D_S); %Right after interface
            
            for i = 1 %y direction BC zero flux (dcdx=0)
                dr = r(2)-r(1);
                dYdt(i,:) = 6*D_C.*(Y(i+1,:)-Y(i,:))./(dr.^2);
                
                j(i,i) = -6*D_C / (dr.^2);
                j(i,i+1) = 6*D_C / (dr.^2);
            end
            for i = 2:ind_Y_r_cavity-1 %inside cavity
                dr = r(i+1)-r(i);
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
                %dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(2*dr); %central difference
                %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
                
                j(i,i-1) = -2*D_C / (r(i)*dr) + D_C / (dr.^2) + v(i) / dr;
                j(i,i) = 2*D_C / (r(i)*dr) - 2*D_C / (dr.^2) - (2/r(i))*v(i) - dvdr(i) - v(i) / dr;
                j(i,i+1) = D_C / (dr.^2);
            end
            for i = ind_Y_r_cavity %right before interface - cavity/stroma
                dr = r(i+1)-r(i);
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
                %dYdr(i,:) = (C_intf_CSa-Y(i-1,:))./(2*dr); %central difference
                %dYdr(i,:) = (C_intf_CSa-Y(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i-1,:)-2.*Y(i,:)+C_intf_CSa)./(dr.^2);
                dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
                
                j(i,i-1) = -2*D_C / (r(i)*dr) + D_C / (dr.^2) + v(i) / dr;
                j(i,i) = 2*D_C / (r(i)*dr) - 2*D_C / (dr.^2) + ...
                    ((D_C^2)/(dr^2))*(1 / ((D_S/P) + D_C)) - ...
                    (2/r(i))*v(i) - dvdr(i) - v(i) / dr;
                j(i,i+1) = ((D_C)/(dr^2))*(D_S / ((D_S/P) + D_C));
            end
            for i = (ind_Y_r_cavity+1) %right after interface - cavity/stroma
                dr = r(i+1)-r(i);
                dYdr(i,:) = (Y(i,:)-C_intf_CSb)./dr; %backward difference
                %dYdr(i,:) = (Y(i+1,:)-C_intf_CSb)./(2*dr); %central difference
                %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr; %forward diffe~rence
                d2Ydr2(i,:) = ((C_intf_CSb)-2.*Y(i,:)+Y(i+1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
                
                j(i,i-1) = 0;
                j(i,i) = (2*D_S / (r(i)*dr)) * (1 - (D_C / (D_C*P + D_S))) + ...
                    ((D_S)/(dr^2))*((D_C / ((D_C*P) + D_S)) - 2) - ...
                    (2/r(i))*v(i) - dvdr(i) - ...
                    (v(i) / dr) * (1 - (D_C / ((D_C*P) + D_S)));
                j(i,i+1) = -(2*(D_S^2) / (r(i))) * (1 / (D_C*P + D_S)) + ...
                    (D_S/(dr^2))*(1 - (D_S / (D_C*P + D_S))) + ...
                    (v(i)/dr)*(D_S / (D_C*P + D_S));
            end
            for i = ind_Y_r_cavity + 2 : (numr - 1) %in stroma
                dr = r(i+1)-r(i);
                %dYdr(i) = (Y(i)-Y(i-1))./dr; %backward difference
                dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(2*dr); %central difference
                %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
                
                j(i,i-1) = -(D_S / (r(i)*dr)) + D_S / (dr^2) + v(i) / (2*dr);
                j(i,i) = -(2*D_S / (dr^2)) - (2/r(i))*v(i) - dvdr(i);
                j(i,i+1) = (D_S / (r(i)*dr)) + D_S / (dr^2) - v(i) / (2*dr);
            
            end
            for i = numr %end of stroma, B.C.
                dr = r(i)-r(i-1);
                %dYdt(i,:) = 6*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2);
                %dYdt(i,:) = D_S*(Y(i,:)-2.*Y(i-1,:)+Y(i-2,:))./(dr.^2)-(2/a)*v(i)*Y(i);
                %dYdt(i,:) = 2*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2);
                dYdt(i,:) = 2*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2) - (2/Rtot)*v(i)*Y(i) - Y(i)*dvdr(i);
                
                j(i-1,i) = (2*D_S / (dr^2));
                j(i,i) = -(2*D_S / (dr^2)) - (2/Rtot)*v(i) - dvdr(i);
            end
            for i = numr+1
                %C_avg_cavity = mean(Y(2:ind_Y_r_cavity,:),1);
                %dYdt(i,:) = (Q_needle * (1 - (theta .* C_oh_inj))) ./ (4 * pi * r0^3 * Y(i,:).^2);
                
                %sum_mass_cavity = 4 * pi * trapz(r(1:ind_Y_r_cavity), (r(1:ind_Y_r_cavity)'.^2) .* Y(1:ind_Y_r_cavity,:), 1);
                %vol_cavity = (4 / 3) * pi * r(ind_Y_r_cavity).^3;
                
                %C_avg_cavity = sum_mass_cavity ./ vol_cavity;
                %dYdt(i,:) = Q_needle * (C_oh_inj - theta .* Y(ind_Y_r_cavity+1,:)) ./ (4 * pi * r0^3 * (kappa + 1) * Y(i,:).^2 .* C_avg_cavity);
                %dYdt(i,:) = Q_needle * (C_oh_inj - epsilon_t .* Y(ind_Y_r_cavity+1,:)) ./ (4 * pi * r0^3 * (kappa + 1) * Y(i,:).^2 .* C_avg_cavity * epsilon_c);
                
                %dYdt(i,:) = strain_helper.concentration_strain_rate(r, Y(1:numr,:), ind_Y_r_cavity, Y(i,:), r0, kappa, Q_needle, C_oh_inj, theta);
                %strain_helper.concentration_strain_rate_porosities_jacobian
                [dlambdadt,jlambda] = helpers{4}(r, Y(1:numr,:), C_intf_CSb, ind_Y_r_cavity, Y(i,:), r_needle, kappa, Q_needle, C_oh_inj, epsilon_c, epsilon_t);
                dYdt(i,:) = dlambdadt;
                j(end,end) = jlambda(end);
                
            end

            F = dYdt;
            J = j;
        end
        function [DYDT] = velocity_bc_cavity_radius_concentration_params(t, Y, r, params, helpers)
            % HERE WE SOLVE FOR THE STRAIN LAMBDA AND THE CONCENTRATION
            % SIMULTANEOUSLY
            %dr = r(2)-r(1);
            dYdt = zeros(numel(r)+1,size(Y,2));
            dYdr = zeros(numel(r)+1,size(Y,2));
            d2Ydr2 = zeros(numel(r)+1,size(Y,2));
%             if exist('ind_r_cavity','var') == 0
%                 ind_r_cavity = ff_helper.binarySearchBin(r, r0);
%             end
            %binarySearchBin = helpers{1};
            %getVelocity = helpers{2};
            %getVelocityGradient = helpers{3};
            %getStrainRate = helpers{4};

            final_ind_r_cavity = params("final_ind_r_cavity");
            r_needle = params("r_needle");
            %theta = params("theta");
            Q_needle = params("Q_needle");
            numr = params("numr");
            D_C = params("D_C");
            D_S = params("D_S");
            P = params("P");
            Rtot = params("Rtot");
            kappa = params("kappa");
            epsilon_t = params("epsilon_t");
            epsilon_c = params("epsilon_c");
            C_oh_inj = params("C_oh_inj");

            if final_ind_r_cavity == 1

                ind_r_cavity = helpers{1}(r, Y(end,end) * r_needle);
                %ind_r_cavity = binarySearchBin(r, Y(end,end) * r_needle);

                if (ind_r_cavity == -1)
                    ind_r_cavity = final_ind_r_cavity;
                end
       
            else
                ind_r_cavity = final_ind_r_cavity;
                
            end
            
            ind_Y_r_cavity = ind_r_cavity;

            %v = getVelocity(r, ind_r_cavity, r_needle, theta, Q_needle);
            %dvdr = getVelocityGradient(r, ind_r_cavity, theta, Q_needle);
            v = helpers{2}(r, ind_r_cavity, r_needle, 1, Q_needle);
            dvdr = helpers{3}(r, ind_r_cavity, 1, Q_needle);
            %v = getVelocity(r, ind_r_cavity, r_needle, 1, Q_needle);
            %dvdr = getVelocityGradient(r, ind_r_cavity, 1, Q_needle);

            %Concentration at cavity/stroma interface
            %C_intf_CSa = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_S.*P)+D_C); %Right before interface
            %C_intf_CSb = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_C./P)+D_S); %Right after interface
            C_intf_CSa = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_S./P)+D_C); %Right before interface
            C_intf_CSb = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_C.*P)+D_S); %Right after interface
            %C_intf_CSa = (D_C.*Y(ind_Y_r_cavity)-(D_C.*Y(ind_Y_r_cavity-1)))./((D_S./P)-D_S); %Right before interface
            %C_intf_CSb = (D_C.*Y(ind_Y_r_cavity)-(D_S.*Y(ind_Y_r_cavity-1)))./((D_S.*P)-D_S); %Right after interface
            
            for i = 1 %y direction BC zero flux (dcdx=0)
                dr = r(2)-r(1);
                dYdt(i,:) = 6*D_C.*(Y(i+1,:)-Y(i,:))./(dr.^2);
            end
            for i = 2:ind_Y_r_cavity-1 %inside cavity
                dr = r(i+1)-r(i);
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
                %dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(2*dr); %central difference
                %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
                %dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
                dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = ind_Y_r_cavity %right before interface - cavity/stroma
                dr = r(i+1)-r(i);
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
                %dYdr(i,:) = (C_intf_CSa-Y(i-1,:))./(2*dr); %central difference
                %dYdr(i,:) = (C_intf_CSa-Y(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i-1,:)-2.*Y(i,:)+C_intf_CSa)./(dr.^2);
                %dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
                dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = (ind_Y_r_cavity+1) %right after interface - cavity/stroma
                dr = r(i+1)-r(i);
                dYdr(i,:) = (Y(i,:)-C_intf_CSb)./dr; %backward difference
                %dYdr(i,:) = (Y(i+1,:)-C_intf_CSb)./(2*dr); %central difference
                %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr; %forward difference
                d2Ydr2(i,:) = ((C_intf_CSb)-2.*Y(i,:)+Y(i+1,:))./(dr.^2);
                %dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = ind_Y_r_cavity + 2 : (numr - 1) %in stroma
                dr = r(i+1)-r(i);
                %dYdr(i) = (Y(i)-Y(i-1))./dr; %backward difference
                dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(2*dr); %central difference
                %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
                %dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = numr %end of stroma, B.C.
                dr = r(i)-r(i-1);
                %dYdt(i,:) = 6*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2);
                %dYdt(i,:) = D_S*(Y(i,:)-2.*Y(i-1,:)+Y(i-2,:))./(dr.^2)-(2/a)*v(i)*Y(i);
                %dYdt(i,:) = 2*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2);
                dYdt(i,:) = 2*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2) - (2/Rtot)*v(i)*Y(i) - Y(i)*dvdr(i);
                %dYdt(i,:) = 0;
            end
            for i = numr+1
                %C_avg_cavity = mean(Y(2:ind_Y_r_cavity,:),1);
                %dYdt(i,:) = (Q_needle * (1 - (theta .* C_oh_inj))) ./ (4 * pi * r0^3 * Y(i,:).^2);
                
                %sum_mass_cavity = 4 * pi * trapz(r(1:ind_Y_r_cavity), (r(1:ind_Y_r_cavity)'.^2) .* Y(1:ind_Y_r_cavity,:), 1);
                %vol_cavity = (4 / 3) * pi * r(ind_Y_r_cavity).^3;
                
                %C_avg_cavity = sum_mass_cavity ./ vol_cavity;
                %dYdt(i,:) = Q_needle * (C_oh_inj - theta .* Y(ind_Y_r_cavity+1,:)) ./ (4 * pi * r0^3 * (kappa + 1) * Y(i,:).^2 .* C_avg_cavity);
                %dYdt(i,:) = Q_needle * (C_oh_inj - epsilon_t .* Y(ind_Y_r_cavity+1,:)) ./ (4 * pi * r0^3 * (kappa + 1) * Y(i,:).^2 .* C_avg_cavity * epsilon_c);
                
                %dYdt(i,:) = strain_helper.concentration_strain_rate(r, Y(1:numr,:), ind_Y_r_cavity, Y(i,:), r0, kappa, Q_needle, C_oh_inj, theta);
                dYdt(i,:) = helpers{4}(r, Y(1:numr,:), dYdr(ind_Y_r_cavity+1,:), ind_Y_r_cavity, Y(i,:), r_needle, kappa, Q_needle, C_oh_inj, epsilon_c, epsilon_t);
                %dYdt(i,:) = getStrainRate(r, Y(1:numr,:), C_intf_CSb, ind_Y_r_cavity, Y(i,:), r_needle, kappa, Q_needle, C_oh_inj, epsilon_c, epsilon_t);

            end

            DYDT = dYdt;
        end
        function [DYDT] = velocity_bc_cavity_radius_concentration(t, Y, r, D_C, D_S, P, theta, epsilon_c, epsilon_t, r_needle, a, numr, Q_needle, kappa, C_oh_inj, ff_helper, v_helper, dvdr_helper, strain_helper, final_ind_r_cavity)
            % HERE WE SOLVE FOR THE STRAIN LAMBDA AND THE CONCENTRATION
            % SIMULTANEOUSLY
            %dr = r(2)-r(1);
            dYdt = zeros(numel(r)+1,size(Y,2));
            dYdr = zeros(numel(r)+1,size(Y,2));
            d2Ydr2 = zeros(numel(r)+1,size(Y,2));
%             if exist('ind_r_cavity','var') == 0
%                 ind_r_cavity = ff_helper.binarySearchBin(r, r0);
%             end

            if final_ind_r_cavity == 1

                ind_r_cavity = ff_helper.binarySearchBin(r, Y(end,end) * r_needle);

                if (ind_r_cavity == -1)
                    ind_r_cavity = final_ind_r_cavity;
                end
       
            else
                ind_r_cavity = final_ind_r_cavity;
                
            end
            
            ind_Y_r_cavity = ind_r_cavity;

            %v = v_helper.ethanol_velocity_porosity(r, ind_r_cavity, r0, theta, Q_needle);
            %dvdr = dvdr_helper.simple_velocity_gradient_porosity(r, ind_r_cavity, theta, Q_needle);
            v = v_helper.ethanol_velocity_porosity(r, ind_r_cavity, r_needle, 1, Q_needle);
            dvdr = dvdr_helper.simple_velocity_gradient_porosity(r, ind_r_cavity, 1, Q_needle);

            %Concentration at cavity/stroma interface
            %C_intf_CSa = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_S.*P)+D_C); %Right before interface
            %C_intf_CSb = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_C./P)+D_S); %Right after interface
            C_intf_CSa = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_S./P)+D_C); %Right before interface
            C_intf_CSb = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_C.*P)+D_S); %Right after interface
            %C_intf_CSa = (D_C.*Y(ind_Y_r_cavity)-(D_C.*Y(ind_Y_r_cavity-1)))./((D_S./P)-D_S); %Right before interface
            %C_intf_CSb = (D_C.*Y(ind_Y_r_cavity)-(D_S.*Y(ind_Y_r_cavity-1)))./((D_S.*P)-D_S); %Right after interface

            
            for i = 1 %y direction BC zero flux (dcdx=0)
                dr = r(2)-r(1);
                dYdt(i,:) = 6*D_C.*(Y(i+1,:)-Y(i,:))./(dr.^2);
            end
            for i = 2:ind_Y_r_cavity-1 %inside cavity
                dr = r(i+1)-r(i);
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
                %dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(2*dr); %central difference
                %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
                %dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
                dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = ind_Y_r_cavity %right before interface - cavity/stroma
                dr = r(i+1)-r(i);
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
                %dYdr(i,:) = (C_intf_CSa-Y(i-1,:))./(2*dr); %central difference
                %dYdr(i,:) = (C_intf_CSa-Y(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i-1,:)-2.*Y(i,:)+C_intf_CSa)./(dr.^2);
                %dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
                dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = (ind_Y_r_cavity+1) %right after interface - cavity/stroma
                dr = r(i+1)-r(i);
                dYdr(i,:) = (Y(i,:)-C_intf_CSb)./dr; %backward difference
                %dYdr(i,:) = (Y(i+1,:)-C_intf_CSb)./(2*dr); %central difference
                %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr; %forward difference
                d2Ydr2(i,:) = ((C_intf_CSb)-2.*Y(i,:)+Y(i+1,:))./(dr.^2);
                %dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = ind_Y_r_cavity + 2 : (numr - 1) %in stroma
                dr = r(i+1)-r(i);
                %dYdr(i) = (Y(i)-Y(i-1))./dr; %backward difference
                dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(2*dr); %central difference
                %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
                %dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = numr %end of stroma, B.C.
                dr = r(i)-r(i-1);
                %dYdt(i,:) = 6*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2);
                %dYdt(i,:) = D_S*(Y(i,:)-2.*Y(i-1,:)+Y(i-2,:))./(dr.^2)-(2/a)*v(i)*Y(i);
                %dYdt(i,:) = 2*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2);
                dYdt(i,:) = 2*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2) - (2/a)*v(i)*Y(i) - Y(i)*dvdr(i);
                %dYdt(i,:) = 0;
            end
            for i = numr+1
                %C_avg_cavity = mean(Y(2:ind_Y_r_cavity,:),1);
                %dYdt(i,:) = (Q_needle * (1 - (theta .* C_oh_inj))) ./ (4 * pi * r0^3 * Y(i,:).^2);
                
                %sum_mass_cavity = 4 * pi * trapz(r(1:ind_Y_r_cavity), (r(1:ind_Y_r_cavity)'.^2) .* Y(1:ind_Y_r_cavity,:), 1);
                %vol_cavity = (4 / 3) * pi * r(ind_Y_r_cavity).^3;
                
                %C_avg_cavity = sum_mass_cavity ./ vol_cavity;
                %dYdt(i,:) = Q_needle * (C_oh_inj - theta .* Y(ind_Y_r_cavity+1,:)) ./ (4 * pi * r0^3 * (kappa + 1) * Y(i,:).^2 .* C_avg_cavity);
                %dYdt(i,:) = Q_needle * (C_oh_inj - epsilon_t .* Y(ind_Y_r_cavity+1,:)) ./ (4 * pi * r0^3 * (kappa + 1) * Y(i,:).^2 .* C_avg_cavity * epsilon_c);
                
                %dYdt(i,:) = strain_helper.concentration_strain_rate(r, Y(1:numr,:), ind_Y_r_cavity, Y(i,:), r0, kappa, Q_needle, C_oh_inj, theta);
                dYdt(i,:) = strain_helper.concentration_strain_rate_porosities(r, Y(1:numr,:), C_intf_CSb, ind_Y_r_cavity, Y(i,:), r_needle, kappa, Q_needle, C_oh_inj, epsilon_c, epsilon_t);
            end

            DYDT = dYdt;
        end
        function [DYDT] = pressure_bc(t, Y, r, D_C, D_S, E, P, theta, k_t, r0, r_tot, numr, Q_needle, C_oh_inj, ff_helper, v_helper, dvdr_helper, p_helper, final_ind_r_cavity)
            % v1.0
            % HERE WE SOLVE FOR THE STRAIN LAMBDA AND THE CONCENTRATION
            % SIMULTANEOUSLY, USING A CONTINUOUS PRESSURE B.C. @ R = R_C
            dr = r(2)-r(1);
            dYdt = zeros(length(r)+1,size(Y,2));
            dYdr = zeros(length(r)+1,size(Y,2));
            d2Ydr2 = zeros(length(r)+1,size(Y,2));

            if final_ind_r_cavity == 1

                ind_r_cavity = ff_helper.binarySearchBin(r, Y(end,end) * r0);

                if (ind_r_cavity == -1)
                    ind_r_cavity = final_ind_r_cavity;
                end
       
            else
                ind_r_cavity = final_ind_r_cavity;
                
            end
            
            ind_Y_r_cavity = ind_r_cavity;

            p_c = p_helper.p_c(E, Y(end,end));
            v = v_helper.pressure_bc(r, ind_r_cavity, r0, r_tot, k_t, p_c, Q_needle, Y(end, end));
            dvdr = dvdr_helper.pressure_bc(r, p_c, r0, r_tot, k_t, Q_needle, Y(end, end), ind_r_cavity);
            
            %Concentration at cavity/stroma interface
            C_intf_CSa = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_S.*P)+D_C); %Right before interface
            C_intf_CSb = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_C./P)+D_S); %Right after interface
            
            for i = 1 %y direction BC zero flux (dcdx=0)
                dYdt(i,:) = 6*D_C.*(Y(i+1,:)-Y(i,:))./(dr.^2);
            end
            for i = 2:ind_Y_r_cavity-1 %inside cavity
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
                %dCdr(i,:) = (C(i+1,:)-C(i-1,:))./(2*dr); %central difference
                %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = ind_Y_r_cavity %right before interface - cavity/stroma
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
                %dCdr(i,:) = (C_intf_CSa-C(i-1,:))./(2*dr); %central difference
                %dCdr(i,:) = (C_intf_CSa-C(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i-1,:)-2.*Y(i,:)+C_intf_CSa)./(dr.^2);
                dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = (ind_Y_r_cavity+1) %right after interface - cavity/stroma
                dYdr(i,:) = (Y(i,:)-C_intf_CSb)./dr; %backward difference
                %dCdr(i,:) = (C(i+1,:)-C_intf_CSb)./(2*dr); %central difference
                %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward diffe~rence
                d2Ydr2(i,:) = ((C_intf_CSb)-2.*Y(i,:)+Y(i+1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = ind_Y_r_cavity+2:(numr - 1) %in stroma
                %dCdr(i) = (C(i)-C(i-1))./dr; %backward difference
                dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(2*dr); %central difference
                %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = numr %end of stroma, B.C.
                %dYdt(i,:) = 6*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2);
                %dYdt(i,:) = D_S*(Y(i,:)-2.*Y(i-1,:)+Y(i-2,:))./(dr.^2)-(2/r_tot)*v(i)*Y(i);

                dYdt(i,:) = 2*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2);
            end
            for i = numr+1
                %sum_mass_cavity = 4 * pi * trapz(r(1:ind_Y_r_cavity), (r(1:ind_Y_r_cavity)'.^2) .* Y(1:ind_Y_r_cavity,:), 1);
                %vol_cavity = (4 / 3) * pi * r(ind_Y_r_cavity).^2;
                %C_avg_cavity = sum_mass_cavity ./ vol_cavity;
                %dYdt(i,:) = Q_needle * (C_oh_inj - theta .* Y(ind_Y_r_cavity+1,:)) ./ (4 * pi * r0^3 * (kappa + 1) * Y(i,:).^2 .* C_avg_cavity);
                %dYdt(i,:) = (Q_needle * (1 - (theta .* C_oh_inj))) ./ (4 * pi * r0^3 * Y(i,:).^2);
                dYdt(i,:) = (Q_needle ./ (4*pi*r0^3*Y(i,:).^2)) - ((k_t * p_c * C_oh_inj) ./ (r0^3 * (Y(i,:).^2))) .* (r_tot ./ (r_tot ./ (Y(i,:) * r0) - 1));
            end

            DYDT = dYdt;
         end
         function [DYDT, ind_r_cavity] = velocity_bc_spline_porosity(t, Y, r, D_C, D_S, P, theta, r0, a, numr, Q_needle, cavity_radius_spline, ff_helper, v_helper, dvdr_helper)
            % HERE WE SOLVE FOR THE STRAIN LAMBDA AND THE CONCENTRATION
            % SIMULTANEOUSLY
            dr = r(2)-r(1);
            dYdt = zeros(length(r),size(Y,2));
            dYdr = zeros(length(r),size(Y,2));
            d2Ydr2 = zeros(length(r),size(Y,2));
                        
            cavity_radius = ppval(cavity_radius_spline,t);
            
            ind_r_cavity = ff_helper.binarySearchBin(r,cavity_radius);
            if (ind_r_cavity == -1)
                ind_r_cavity = numel(r);
            end
          
            v = v_helper.ethanol_velocity_porosity(r, ind_r_cavity, r0, theta, Q_needle);
            dvdr = dvdr_helper.simple_velocity_gradient_porosity(r, ind_r_cavity, theta, Q_needle);
            
            %Concentration at cavity/stroma interface
            C_intf_CSa = (D_C.*Y(ind_r_cavity)+(D_S.*Y(ind_r_cavity+1)))./((D_S.*P)+D_C); %Right before interface
            C_intf_CSb = (D_C.*Y(ind_r_cavity)+(D_S.*Y(ind_r_cavity+1)))./((D_C./P)+D_S); %Right after interface
            
            for i = 1 %y direction BC zero flux (dcdx=0)
                dYdt(i,:) = 6*D_C.*(Y(i+1,:)-Y(i,:))./(dr.^2);
            end
            for i = 2:ind_r_cavity-1 %inside cavity
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
                %dCdr(i,:) = (C(i+1,:)-C(i-1,:))./(2*dr); %central difference
                %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = ind_r_cavity %right before interface - cavity/stroma
               dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
               %dCdr(i,:) = (C_intf_CSa-C(i-1,:))./(2*dr); %central difference
               %dCdr(i,:) = (C_intf_CSa-C(i,:))./dr; %forward difference
               d2Ydr2(i,:) = (Y(i-1,:)-2.*Y(i,:)+C_intf_CSa)./(dr.^2);
               dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = (ind_r_cavity+1) %right after interface - cavity/stroma
                dYdr(i,:) = (Y(i,:)-C_intf_CSb)./dr; %backward difference
                %dCdr(i,:) = (C(i+1,:)-C_intf_CSb)./(2*dr); %central difference
                %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward diffe~rence
                d2Ydr2(i,:) = ((C_intf_CSb)-2.*Y(i,:)+Y(i+1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = ind_r_cavity+2:(numr-1) %in stroma
                %dCdr(i) = (C(i)-C(i-1))./dr; %backward difference
                dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(2*dr); %central difference
                %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = numr %end of stroma, B.C.
                %dCdt(i,:) = 6*D_S*(C(i-1,:)-C(i,:))./(dr.^2);
                dYdt(i,:) = D_S*(Y(i,:)-2.*Y(i-1,:)+Y(i-2,:))./(dr.^2)-(2/a)*v(i)*Y(i);

                %dCdt(i,:) = 2*D_S*(C(i-1,:)-C(i,:))./(dr.^2);
            end

            DYDT = dYdt;
        end
        function [DYDT] = velocity_bc_spline(t, Y, r, D_C, D_S, phi, a, numr, v, dvdr, cavity_radius_spline, helper)
            % HERE WE SOLVE FOR THE STRAIN LAMBDA AND THE CONCENTRATION
            % SIMULTANEOUSLY
            dr = r(2)-r(1);
            dYdt = zeros(length(r),size(Y,2));
            dYdr = zeros(length(r),size(Y,2));
            d2Ydr2 = zeros(length(r),size(Y,2));
            
            % cavity_radius_t = spline(tvec, cavity_radius, t);
            
            cavity_radius = ppval(cavity_radius_spline,t);
            
            ind_r_cavity = helper.binarySearchBin(r,cavity_radius);
            if (ind_r_cavity == -1)
                ind_r_cavity = numel(r);
            end     

          
%             for i = 1 % here we solve for the strain
%                 dYdt(i,:) = Q_needle*(1-omega)./(4*pi*(r0^3)*(Y(i,:).^2));
%             end
%            
%             cavity_radius = dYdt(1,:) * r0;
% 
%             ind_r_cavity = helper.binarySearchBin(r,cavity_radius(end));
%             if (ind_r_cavity == -1)
%                 ind_r_cavity = numel(r);
%             end
%             ind_Y_r_cavity = ind_r_cavity+1;
            
            %Concentration at cavity/stroma interface
            C_intf_CSa = (D_C.*Y(ind_r_cavity)+(D_S.*Y(ind_r_cavity+1)))./((D_S.*phi)+D_C); %Right before interface
            C_intf_CSb = (D_C.*Y(ind_r_cavity)+(D_S.*Y(ind_r_cavity+1)))./((D_C./phi)+D_S); %Right after interface
            
            for i = 1 %y direction BC zero flux (dcdx=0)
                dYdt(i,:) = 6*D_C.*(Y(i+1,:)-Y(i,:))./(dr.^2);
            end
            for i = 2:ind_r_cavity-1 %inside cavity
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
                %dCdr(i,:) = (C(i+1,:)-C(i-1,:))./(2*dr); %central difference
                %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = ind_r_cavity %right before interface - cavity/stroma
               dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
               %dCdr(i,:) = (C_intf_CSa-C(i-1,:))./(2*dr); %central difference
               %dCdr(i,:) = (C_intf_CSa-C(i,:))./dr; %forward difference
               d2Ydr2(i,:) = (Y(i-1,:)-2.*Y(i,:)+C_intf_CSa)./(dr.^2);
               dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = (ind_r_cavity+1) %right after interface - cavity/stroma
                dYdr(i,:) = (Y(i,:)-C_intf_CSb)./dr; %backward difference
                %dCdr(i,:) = (C(i+1,:)-C_intf_CSb)./(2*dr); %central difference
                %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward diffe~rence
                d2Ydr2(i,:) = ((C_intf_CSb)-2.*Y(i,:)+Y(i+1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = ind_r_cavity+2:(numr-1) %in stroma
                %dCdr(i) = (C(i)-C(i-1))./dr; %backward difference
                dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(2*dr); %central difference
                %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = numr %end of stroma, B.C.
                %dCdt(i,:) = 6*D_S*(C(i-1,:)-C(i,:))./(dr.^2);
                dYdt(i,:) = D_S*(Y(i,:)-2.*Y(i-1,:)+Y(i-2,:))./(dr.^2)-(2/a)*v(i)*Y(i);

                %dCdt(i,:) = 2*D_S*(C(i-1,:)-C(i,:))./(dr.^2);
            end

            DYDT = dYdt;
        end
        function [DYDT] = velocity_bc(t, Y, r, omega, Q_needle, r0, D_C, D_S, phi, a, numr, v, dvdr, helper)
            % HERE WE SOLVE FOR THE STRAIN LAMBDA AND THE CONCENTRATION
            % SIMULTANEOUSLY
            dr = r(2)-r(1);
            dYdt = zeros(length(r)+1,size(Y,2));
            dYdr = zeros(length(r)+1,size(Y,2));
            d2Ydr2 = zeros(length(r)+1,size(Y,2));
            
          
            for i = 1 % here we solve for the strain
                dYdt(i,:) = Q_needle*(1-omega)./(4*pi*(r0^3)*(Y(i,:).^2));
            end
           
            cavity_radius = dYdt(1,:) * r0;

            ind_r_cavity = helper.binarySearchBin(r,cavity_radius(end));
            if (ind_r_cavity == -1)
                ind_r_cavity = numel(r);
            end
            ind_Y_r_cavity = ind_r_cavity+1;
            
            %Concentration at cavity/stroma interface
            C_intf_CSa = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_S.*phi)+D_C); %Right before interface
            C_intf_CSb = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_C./phi)+D_S); %Right after interface
            
            for i = 2 %y direction BC zero flux (dcdx=0)
                dYdt(i,:) = 6*D_C.*(Y(i+1,:)-Y(i,:))./(dr.^2);
            end
            for i = 3:ind_Y_r_cavity-1 %inside cavity
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
                %dCdr(i,:) = (C(i+1,:)-C(i-1,:))./(2*dr); %central difference
                %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = ind_Y_r_cavity %right before interface - cavity/stroma
               dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
               %dCdr(i,:) = (C_intf_CSa-C(i-1,:))./(2*dr); %central difference
               %dCdr(i,:) = (C_intf_CSa-C(i,:))./dr; %forward difference
               d2Ydr2(i,:) = (Y(i-1,:)-2.*Y(i,:)+C_intf_CSa)./(dr.^2);
               dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = (ind_Y_r_cavity+1) %right after interface - cavity/stroma
                dYdr(i,:) = (Y(i,:)-C_intf_CSb)./dr; %backward difference
                %dCdr(i,:) = (C(i+1,:)-C_intf_CSb)./(2*dr); %central difference
                %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward difference
                d2Ydr2(i,:) = ((C_intf_CSb)-2.*Y(i,:)+Y(i+1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = ind_Y_r_cavity+2:(numr-1) %in stroma
                %dCdr(i) = (C(i)-C(i-1))./dr; %backward difference
                dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(2*dr); %central difference
                %dCdr(i,:) = (C(i+1,:)-C(i,:))./dr; %forward difference
                d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v(i).*Y(i,:)))-(Y(i,:).*dvdr(i))-(v(i).*dYdr(i,:));
            end
            for i = numr %end of stroma, B.C.
                %dCdt(i,:) = 6*D_S*(C(i-1,:)-C(i,:))./(dr.^2);
                dYdt(i,:) = D_S*(Y(i,:)-2.*Y(i-1,:)+Y(i-2,:))./(dr.^2)-(2/a)*v(i)*Y(i);

                %dCdt(i,:) = 2*D_S*(C(i-1,:)-C(i,:))./(dr.^2);
            end

            DYDT = dYdt;
        end
    end
end
