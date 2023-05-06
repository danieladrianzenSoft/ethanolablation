classdef fastFindLib
    methods(Static)
        function [idx] = binarySearch(sorted_array, num)
            left = 1;
            right = length(sorted_array);
            flag = 0;

            while (left <= right)
                mid = ceil((left + right) / 2);

                if (sorted_array(mid) == num)
                    idx = mid;
                    flag = 1;
                    break;
                elseif (sorted_array(mid) > num)
                    right = mid - 1;
                else
                    left = mid + 1;
                end
            end

            if (flag == 0)
                idx = -1;
            end
        end

        function [idx] = binarySearchBin(sorted_array, num)
            left = 1;
            right = length(sorted_array);
            flag = 0;

            while (right - left > 1)
                mid = ceil((left + right) / 2);

                if (num == sorted_array(mid))
                    idx = mid;
                    flag = 1;
                    break;
                elseif (num < sorted_array(mid))
                    right = mid;
                else
                    left = mid;
                end
            end

            if (flag == 0)
                if (sorted_array(left) <= num && num < sorted_array(right))
                    idx = left;
                else
                    idx = -1;
                end
            end
        end

        function [varargout] = fastIntersect(a, b, option)
        % FASTINTERSECT    Optimized intersection on sorted inputs.
        %
        %    This function performs the intersection operation on two SORTED
        %    one-dimensional arrays. It is important to note that this function
        %    will respect non-unique values in either of the input arrays. Hence,
        %    it does not perform a strict set intersection. Furthermore, this
        %    function requires a third argument, one that specifies what this
        %    function will return. The options are as follows:
        %
        %    'size'      - Return the size of the intersection.
        %    'bool'      - Return true if the intersection is nonempty.
        %    'bool vec'  - Return booleon vector (size of first argument) that
        %                  labels elements in the intersection (including repeats);
        %                  including 2nd output will swap 1st & 2nd argument order.
        %    'all'       - Return true if all elements of the smaller array are
        %                  contained within the larger array.
        %    'elements'  - Return a one-dimensional array containing the elements
        %                  in the intersection.
        %    'indices'   - Return the indices of the elements in the first argument
        %                  that are in the intersection (including repeats);
        %                  including 2nd output will swap 1st & 2nd argument order.
        %    'index'     - Same as the above.
        %
        % Author: Petes

            if strcmp(option, 'size')
                varargout{1} = fastFindLib.fastIntersectSize(a, b);
            elseif strcmp(option, 'bool')
                varargout{1} = fastFindLib.fastIntersectBool(a, b);
            elseif strcmp(option, 'bool vec') || strcmp(option, 'boolvec') || strcmp(option, 'bool vector')
                varargout{1} = fastFindLib.fastIntersectBoolVec(a, b);
                if nargout == 2
                    varargout{2} = fastFindLib.fastIntersectBoolVec(b, a);
                end
            elseif strcmp(option, 'all')
                varargout{1} = fastFindLib.fastIntersectAll(a, b);
            elseif strcmp(option, 'elements')
                varargout{1} = fastFindLib.fastIntersectElements(a, b);
            elseif strcmp(option, 'indices') || strcmp(option, 'index')
                varargout{1} = fastFindLib.fastIntersectIndices(a, b);
                if nargout == 2
                    varargout{2} = fastFindLib.fastIntersectIndices(b, a);
                end
            else
                error('Error: Option string unrecognized.');
            end
        end

        function [sz] = fastIntersectSize(a, b)
            aEnd = length(a);
            bEnd = length(b);
            aIter = 1;
            bIter = 1;
            sz = 0;
            while (aIter <= aEnd && bIter <= bEnd)
                if (a(aIter) < b(bIter))
                    aIter = aIter + 1;
                elseif (a(aIter) > b(bIter))
                    bIter = bIter + 1;
                else
                    aIter = aIter + 1;
                    bIter = bIter + 1;
                    sz = sz + 1;
                end
            end
        end

        function [flag] = fastIntersectBool(a, b)
            aEnd = length(a);
            bEnd = length(b);
            aIter = 1;
            bIter = 1;
            flag = false;
            while (aIter <= aEnd && bIter <= bEnd)
                if (a(aIter) < b(bIter))
                    aIter = aIter + 1;
                elseif (a(aIter) > b(bIter))
                    bIter = bIter + 1;
                else
                    flag = true;
                    break;
                end
            end
        end

        function [idxA] = fastIntersectBoolVec(a, b)
            aEnd = length(a);
            bEnd = length(b);
            idxA = false(aEnd, 1);
            aIter = 1;
            bIter = 1;
            while (aIter <= aEnd && bIter <= bEnd)
                if (a(aIter) < b(bIter))
                    aIter = aIter + 1;
                elseif (a(aIter) > b(bIter))
                    bIter = bIter + 1;
                else
                    currValue = a(aIter);
                    while (aIter <= aEnd && a(aIter) == currValue)
                        idxA(aIter) = true;
                        aIter = aIter + 1;
                    end
                    while (bIter <= bEnd && b(bIter) == currValue)
                        bIter = bIter + 1;
                    end
                end
            end
        end

        function [out] = fastIntersectAll(a, b)
            aEnd = length(a);
            bEnd = length(b);
            % force a to be the larger array
            if bEnd > aEnd
                [a, b] = swap(a, b);
                [aEnd, bEnd] = swap(aEnd, bEnd);
            end
            aIter = 1;
            bIter = 1;
            while (aIter <= aEnd && bIter <= bEnd)
                if (a(aIter) < b(bIter))
                    aIter = aIter + 1;
                elseif (a(aIter) > b(bIter))
                    out = false;
                    return;
                else
                    aIter = aIter + 1;
                    bIter = bIter + 1;
                end
            end
            out = (bIter - 1 == bEnd);
        end

        function [out] = fastIntersectElements(a, b)
            aEnd = length(a);
            bEnd = length(b);
            out = zeros(max(aEnd, bEnd), 1);
            aIter = 1;
            bIter = 1;
            oIter = 1;
            while (aIter <= aEnd && bIter <= bEnd)
                if (a(aIter) < b(bIter))
                    aIter = aIter + 1;
                elseif (a(aIter) > b(bIter))
                    bIter = bIter + 1;
                else
                    out(oIter) = a(aIter);
                    aIter = aIter + 1;
                    bIter = bIter + 1;
                    oIter = oIter + 1;
                end
            end
            % Over allocated earlier, trim off the end
            out = out(1 : oIter - 1);
        end

        function [idxA] = fastIntersectIndices(a, b)
            aEnd = length(a);
            bEnd = length(b);
            idxA = zeros(aEnd, 1);
            aIter = 1;
            bIter = 1;
            aOutIter = 1;
            while (aIter <= aEnd && bIter <= bEnd)
                if (a(aIter) < b(bIter))
                    aIter = aIter + 1;
                elseif (a(aIter) > b(bIter))
                    bIter = bIter + 1;
                else
                    currValue = a(aIter);
                    while (aIter <= aEnd && a(aIter) == currValue)
                        idxA(aOutIter) = aIter;
                        aOutIter = aOutIter + 1;
                        aIter = aIter + 1;
                    end
                    while (bIter <= bEnd && b(bIter) == currValue)
                        bIter = bIter + 1;
                    end
                end
            end
            % Output non-zeros
            idxA = idxA(1 : aOutIter - 1);
        end

        function [b, a] = swap(a, b)
        end
    end
end