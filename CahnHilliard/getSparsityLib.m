classdef getSparsityLib
    methods(Static)
        function S = cahn_hilliard_2D(Nx,Ny)
            e = ones(Ny,1);
            A = spdiags([e e e e e],-2:2,Ny,Ny);
            v = full(A);
            S = repmat(v,Nx,Nx);
        end
        function S = cahn_hilliard_1D(Nx,Ny)
            totalSize = Nx;
            S = ones(totalSize,totalSize);
        end
    end
end
