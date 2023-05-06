function [inds] = mask2ind(mask)
    [Ny,Nx] = size(mask);
    [x_vals,y_vals] = find(mask>0);
    %coords = [x_vals,y_vals];
    inds = sub2ind([Ny,Nx],y_vals,x_vals);
end