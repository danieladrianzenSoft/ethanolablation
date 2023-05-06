function neighbours = getNeighbours(perimeter,coords)
    [r,c]=size(perimeter);
    checkMatrix = [0,-1; %down
        0,1; %up
        1,0; %right
        -1,0; %left
        1,1; %right-up
        1,-1; %right-down
        -1,1; %left-up
        -1,-1]; %left-down
    for i=1:size(checkMatrix,1)
        if (coords(1) ~= 1 && coords(1) ~= c && coords(2) ~= 1 && coords(2) ~= r)
            neighbours = repmat(coords,size(checkMatrix,1),1) + checkMatrix;
        elseif coords(1) == 1
            if (coords(2) == 1)
                neighbours = repmat(coords,3,1) + checkMatrix([2,3,5],:);
            elseif (coords(2) == r)
                neighbours = repmat(coords,3,1) + checkMatrix([1,3,6],:);
            else
                neighbours = repmat(coords,5,1) + checkMatrix([1,2,3,5,6],:);
            end
        elseif coords(1) == c
            if (coords(2) == 1)
                neighbours = repmat(coords,3,1) + checkMatrix([2,4,7],:);
            elseif (coords(2) == r)
                neighbours = repmat(coords,3,1) + checkMatrix([1,4,8],:);
            else
                neighbours = repmat(coords,5,1) + checkMatrix([1,2,4,7,8],:);
            end
        elseif coords(2) == 1
            neighbours = repmat(coords,5,1) + checkMatrix([2,3,4,5,7],:);
        elseif coords(2) == r
            neighbours = repmat(coords,5,1) + checkMatrix([1,3,4,6,8],:);
        end
    end
end