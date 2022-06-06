function [T] = getSolutionTransitionMatrix(boutMatrices)
T = zeros(4);
for i=1:4
    for j=i:4
        if (i == 1 && j > 1)
            continue;
        elseif (i == j)
            T(i,j) = mean([boutMatrices{i,j}(2,1) boutMatrices{i,j}(1,2)]);
        else
            T(j,i) = boutMatrices{i,j}(2,1);
            T(i,j) = boutMatrices{i,j}(1,2);
        end
    end
end
end

