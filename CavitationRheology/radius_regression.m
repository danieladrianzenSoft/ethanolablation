function rc_out = radius_regression(r_fit,teval,t,cavity_radius)
    
    current_r = zeros(length(teval),1);
    n_tmax = find(teval >= t(end),1);
    

    if (isempty(n_tmax))
%         if (length(t) == 1)
%            current_r(1) = cavity_radius(end);
%         else
        current_r = polyval(r_fit,teval);
%         end
    else
        current_r(1 : n_tmax) = polyval(r_fit,teval(1:n_tmax));
        current_r(n_tmax + 1 : end) = cavity_radius(end);
    end
    %current_r = polyfit(t,cavity_radius,5);
    rc_out = current_r;
end