classdef getHelperFunctionLib
    methods(Static)
        function v = getVelocity_velocity_bc_porosity(r, ind_r_cavity, r0, phi, Q_needle, getVelocityLib)
            v = getVelocityLib.ethanol_velocity_v2(r, ind_r_cavity, r0, phi, Q_needle);
        end
        function v = getVelocityGradient_velocity_bc_porosity(r, ind_r_cavity, r0, phi, Q_needle, getVelocityLib)
            v = getVelocityLib.ethanol_velocity_v2(r, ind_r_cavity, r0, phi, Q_needle);
        end
    end
end
