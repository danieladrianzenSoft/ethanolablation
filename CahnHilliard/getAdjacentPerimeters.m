function [oPerimeter,iPerimeter] = getAdjacentPerimeters(matrix,perimeter)
    [xinds,yinds] = find(perimeter>0);
    inds = [xinds,yinds];
    oPerimeter = zeros(size(matrix));
    iPerimeter = zeros(size(matrix));
    %oPerimInds = [];
    %iPerimInds = [];
    for i = 1:size(inds,1)
        neighbours = getNeighbours(perimeter,inds(i,:));
        for j = 1:size(neighbours,1)
            neighbour = neighbours(j,:);
            if matrix(neighbour(1),neighbour(2)) > 0
                iPerimeter(neighbour(1),neighbour(2)) = 1;
            elseif matrix(neighbour(1),neighbour(2)) < 0
                oPerimeter(neighbour(1),neighbour(2)) = 1;
            end
        end
    end
end