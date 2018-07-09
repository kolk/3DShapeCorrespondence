function [M, distMatrix, l] = associationGraphSymmetry(vertices1, vertices2, shape1, shape2)
    %vertices : Map of graphs of symmetries in shapes
    %vertices2: graph of symmetries in shape2
    %shape1: mesh for shape1
    %shape2: mesh for shape2
    geoDist1 = findGeodesic(vertices1, shape1);
    geoDist2 = findGeodesic(vertices2, shape2);

    n1 = length(vertices1);
    n2 = length(vertices2);
    n = n1 + n2;
    l = [vertices1; vertices2];
    M = zeros(n,n);
    Infinity = 1000000;
    distMatrix = zeros(n,n);

   % Affinity Matrix has non-0 weights only on the off-diagonal
   % elements
    for i = 1:n1
        for j = n1+1:n
            distMatrix(i,j) = abs(geoDist1(i) - geoDist2(j-n1))
        end
    end
   for i = n1+1:n1+n1
        for j=1:n2
            distMatrix(i,j) = abs(geoDist1(i-n1) - geoDist2(j))
        end
    end



    for i = 1:n
        for j = 1:n
            if(distMatrix(i,j) < Infinity)
                if(distMatrix(i,j) > 0)
                    M(i,j) = exp(-distMatrix(i,j)/(mean(mean(distMatrix))));
                end
            end
        end
    end
end