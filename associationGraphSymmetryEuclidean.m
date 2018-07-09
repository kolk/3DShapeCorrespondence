function [M, distMatrix, l] = associationGraphSymmetryEuclidean(vertices, data, features, interestPoints)
    %vertices : vertices of symmetries 
    %data: data points of shape in embedding space
    % function returns the Association graph matrix M, distance matrix
    % distMatrix and all possible associations l. M is found by considering
    % the euclidean distance of symmetric points in embedding space.

    n = length(vertices);
    l = vertices;
    M = zeros(n,n);
    distMatrix = zeros(n, n);
    distDiag = zeros(n,1);
    pairFeatDist = zeros(n, n);
    
    for i = 1:n
        for j = 1:n
            % compare features 
            ipIndx1 = find(interestPoints == vertices(i,1));
            ipIndx2 = find(interestPoints == vertices(i,2));
            euclideanDist = sqrt(sum(features(ipIndx1,:) .* features(ipIndx2,:)));
            pairFeatDist(i,j) = euclideanDist;
        end
    end
    
   % Affinity Matrix has non-0 weights only on the off-diagonal
   % elements
    for i = 1:n
        for j = i:n
            % for off diagonal elements compute euclidean distance of the
            % symmetric points
            if i ~= j
                euclideanDist1 = sqrt(sum((data(vertices(i,1),:) - data(vertices(j,1),:)).^2))
                euclideanDist2 = sqrt(sum((data(vertices(i,2),:) - data(vertices(j,2),:)).^2))
                dist = abs(euclideanDist1 - euclideanDist2)
                distMatrix(i,j) = dist;
                distMatrix(j,i) = distMatrix(i,j);
            else
                % for diagonal elements compare the features by bhatt dist
                ipIndx1 = find(interestPoints == vertices(i,1));
                ipIndx2 = find(interestPoints == vertices(i,2));
                euclideanDist = sqrt(sum(features(ipIndx1,:) .* features(ipIndx2,:)));
                distDiag(i,j) = euclideanDist;
            end
            

            
        end
    end

    
    
    for i = 1:n
        for j = 1:n
            if i == j
                M(i,i) = exp(-distDiag(i,i)/mean(mean(distDiag)));
            else
               if(distMatrix(i,j) > 0)
                    M(i,j) = exp(-distMatrix(i,j)/(mean(mean(distMatrix))));
               end
            end
        end
    end
end