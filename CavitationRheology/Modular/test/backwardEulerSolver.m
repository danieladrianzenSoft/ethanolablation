function [x,y] = backwardEulerSolver(f, y0, xSpan, dx)
    numX = (xSpan(2)-xSpan(1)) / dx;
    
    x = [xSpan(1), zeros(1,numX)];
    y = [y0, zeros(length(y0), numX)];
    
    for iter = 1:numX
        x(iter+1) = x(iter) + dx;
        y_next = y(:,iter) + dx * f(x(iter), y(:,iter));
        y(:,iter+1) = y(:,iter) + dx * f(x(iter + 1), y_next);
    end
end