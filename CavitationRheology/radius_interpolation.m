function rc_out = radius_interpolation(teval,t,cavity_radius)
    
    current_r = zeros(length(teval),1);
    n_tmax = find(teval >= t(end),1);
    

    if (isempty(n_tmax))
%         if (length(t) == 1)
%            current_r(1) = cavity_radius(end);
%         else
        current_r = interp1(t,cavity_radius,teval,'spline');
%         end
    else
        current_r(1 : n_tmax) = interp1(t,cavity_radius,teval(1 : n_tmax),'spline');
        current_r(n_tmax + 1 : end) = cavity_radius(end);
    end
    %current_r = polyfit(t,cavity_radius,5);
    rc_out = current_r;
end