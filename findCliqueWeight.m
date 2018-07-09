function weight = findCliqueWeight(M, correspondences, numGaussians)
    weight = 0;
    for i = 1:size(correspondences,1)
        for j = i+1:size(correspondences,1)
           index1 = ((correspondences(i,1) - 1) * numGaussians) + correspondences(i,2);  
           index2 = ((correspondences(j,1) - 1) * numGaussians) + correspondences(j,2);  
           weight = weight + M(index1, index2);
        end
    end
end