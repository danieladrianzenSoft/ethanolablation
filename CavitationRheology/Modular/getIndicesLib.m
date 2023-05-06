classdef getIndicesLib
    methods(Static)
        function ind = getIndNeedleRadius(r, r0)
            ff = fastFindLib;
            ind = ff.binarySearchBin(r, r0);
        end
        function ind = getIndCavityRadius(t, r, cavity_radius)
            % get index in r vector of cavity radius over time
            ff = fastFindLib;

            ind = zeros(1,length(t)); 

            for jj = 1:length(ind)
               ind(jj) = ff.binarySearchBin(r,cavity_radius(jj));
            end

        end
        function ind = getIndInjectionEnd(t, vol, q)
            % get index in t of the time when the injection ends 
            
            ff = fastFindLib;

            t_injection_end = vol/q;

            ind = ff.binarySearchBin(t, t_injection_end);

        end
    end
end