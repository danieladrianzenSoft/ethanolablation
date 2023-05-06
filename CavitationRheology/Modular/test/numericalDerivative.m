function df = numericalDerivative(f, x, y, dx, wrt)
    if wrt == 'x'
        df = (f(x + dx,y) - f(x,y)) / dx;
    else
        df = (f(x,y+dx) - f(x,y)) / dx;
    end
end