function [x,y] = midpointImplicitSolver(f, J, y0, xSpan, dx, params, helpers, options)

%     default_tolerance = 1e-6;
%     default_maxIters = 200;
% 
%     p = inputParser;
%     addOptional(p, 'tolerance', default_tolerance);
%     addOptional(p, 'maxIters', default_maxIters); 
%             
%     parse(p,varargin{:});

    if isKey(options,"tolerance") == 0
       options("tolerance") = 1e-6;
    end
    if isKey(options,"maxIters") == 0
       options("maxIters") = 200; 
    end
    
    numX = floor((xSpan(2)-xSpan(1)) / dx) + 1;
    
    x = [xSpan(1), zeros(1,numX)];
    y = [y0, zeros(length(y0), numX)];
    
    for i = 1:numX
        x(i+1) = x(i) + dx;
        %y_guess = y(:,i) + dx * f(x(i), y(:,i));
        y_guess = y(:,i) + dx * f(x(i) + dx/2, y(:,i) + (dx/2) * (f(x(i),y(:,i))));
        y(:,i+1) = newtonRaphson(f, J, params, x(i), x(i+1), y(:,i), y_guess, dx, options("tolerance"), options("maxIters"), helpers);
    end
end

function Y = newtonRaphson(f, J, params, xcurr, xnext, ycurr, Y, dx, tolerance, maxIters, helpers)
    error = 1;
    %maxIters = 200;
    iter = 1;
    while (error >= tolerance && iter <= maxIters)
        F = midpointStep(xcurr, xnext, ycurr, Y, dx, f);
        %params.y = ycurr;
        y_j = (1/2) * (ycurr + Y);
        x_j = (xcurr + dx/2);
        jacobian_next = J(params, y_j, x_j, helpers);
        diff = ((dx)*(jacobian_next)-eye(numel(Y)))\F;
        %diff = (dx)*(jacobian_next-eye(numel(Y)))\F;
        Y = Y - diff;
        %error = abs(y_new - Y)/y_new;
        %error = (norm(y_new, inf) - norm(Y,inf)) / norm(Y, inf);
        error = norm(diff, inf) / norm(Y,inf);
        %Y = y_new;
        iter = iter + 1;
    end
end

function F = midpointStep(xcurr, xnext, ycurr, ynext, dx, f)
    F = ycurr + (dx) * (f(xcurr + dx/2, (1/2) * (ycurr + ynext))) - ynext;
end

