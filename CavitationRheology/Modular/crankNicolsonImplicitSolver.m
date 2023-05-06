function [x,y] = crankNicolsonImplicitSolver(f, y0, xSpan, dx, options)

    if isKey(options,"tolerance") == 0
       options("tolerance") = 1e-6;
    end
    if isKey(options,"maxIters") == 0
       options("maxIters") = 200; 
    end
    
    numX = floor((xSpan(2)-xSpan(1)) / dx) + 1;

    x = [xSpan(1); zeros(numX,1)];
    y = [y0, zeros(length(y0), numX)];
    
    for i = 1:numX
        x(i+1) = x(i) + dx;
        y_guess = y(:,i) + dx * f(x(i), y(:,i));
        y(:,i+1) = newtonRaphson(f, x(i), x(i+1), y(:,i), y_guess, dx, options("tolerance"), options("maxIters"));
    end
end

function Y = newtonRaphson(f, xcurr, xnext, ycurr, Y, dx, tolerance, maxIters)
    error = 1;
    iter = 1;
    while (error >= tolerance && iter <= maxIters)
        [F, jcurr, jnext] = crankNicolsonStep(xcurr, xnext, ycurr, Y, dx, f);
        diff = ((dx/2)*(jcurr + jnext)-eye(numel(Y)))\F;
        Y = Y - diff;
        error = norm(diff, inf) / norm(Y,inf);
        iter = iter + 1;
    end
end

function [F, jcurr, jnext] = crankNicolsonStep(xcurr, xnext, ycurr, ynext, dx, f)
    [fnext, jnext] = f(xnext, ynext);
    [fcurr, jcurr] = f(xcurr, ycurr); 
    F = ycurr + (dx / 2) * (fnext + fcurr) - ynext;
end
