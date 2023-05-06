classdef getSparsityLib
    methods(Static)
         function S = diagonalBandLastElem(numr)
            % diagonal with +1 -1 padding, plus end row
            % and end col, for when calculating
            % lambda simultaneous with concentration.
            totalSize = numr+1;
            e = ones(totalSize,1);
            A = spdiags([e e e],-1:1,totalSize,totalSize);
            S = full(A);
            S(end,end) = 1;
            %S(1:end,end) = 1;
         end
         function S = diagonalBandLastRowLastCol(numr)
            % diagonal with +1 -1 padding, plus end row
            % and end col, for when calculating
            % lambda simultaneous with concentration.
            totalSize = numr+1;
            e = ones(totalSize,1);
            A = spdiags([e e e],-1:1,totalSize,totalSize);
            S = full(A);
            S(end,1:end) = 1;
            S(1:end,end) = 1;
         end
         function S = diagonalBandFirstRowFirstCol(numr)
            % diagonal with +1 -1 padding, plus first row
            % and first col, for when calculating
            % lambda simultaneous with concentration.
            totalSize = numr+1;
            e = ones(totalSize,1);
            A = spdiags([e e e],-1:1,totalSize,totalSize);
            S = full(A);
            S(1,1:end) = 1;
            S(1:end,1) = 1;
        end
        function S = diagonalBand_v2(numr)
            totalSize = numr+1;
            e = ones(totalSize,1);
            A = spdiags([e e e],-1:1,totalSize,totalSize);
            S = full(A);
        end
        function S = diagonalBand(numr)
            totalSize = numr;
            e = ones(totalSize,1);
            A = spdiags([e e e],-1:1,totalSize,totalSize);
            S = full(A);
        end
        function S = all(numr)
            totalSize = numr+1;
            S = ones(totalSize,totalSize);
        end
    end
end




