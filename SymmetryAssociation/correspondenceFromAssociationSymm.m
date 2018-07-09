function [correspondences, indices] = correspondenceFromAssociationSymm(M, l)
    % computer eigenvec and eigenval of M
    [V, D] = eig(M);
    D = diag(D);
    
    % sorted eigen values
    [D, ind] = sort(D,'descend'); 
    
    % sorted eigen vectors
    V = V(:,ind); 
    
    % principal eigenvector
    xstar = V(:,1);
    
    % number of nodes
    n = size(M,1);
    k = n;
    %{
    % find optimum number of symmetries
    for i = 1:n
        for j = 1:n
            if(i ~= j)
                dist = sqrt(sum(features(ipIndx1,:) .* features(ipIndx2,:)));
                feat(i,:) .* feat(j,:))
            end
        end
    end
    %}
    
    
    % solution vector
    sol = zeros(n,1);
    
    % if all values of principal eigenvector -ve, convert to +ve  
    if(sum(find(xstar < 0)) == length(xstar))
        xstar = abs(xstar);
    end
    [sortxstar, indx] = sort(xstar, 'descend');
    
    
    % indices chosen as good correspondes
    indices = zeros(n,1);
   

    invalid = min(xstar) - 0.01;
    disp('calculating correspondences')
    numCorrespondences = 0;
    while numCorrespondences < k
        [~,astar] = max(xstar);

        if (xstar(astar) == invalid)
            %return sol;
            break;
        else
            sol(astar) = 1;
            
            disp([l(astar,:) xstar(astar)])
            %remove from L all potential conflicts with a* = (i, j), ie,
            %assignments of the form (i,k) and (q, j).
            xstar(astar) = invalid;
            i = l(astar, 1);
            j = l(astar, 2);
            iIndx = find(l(:,1) == i);
            jIndx = find(l(:,2) == j);
            xstar(iIndx) = invalid;
            xstar(jIndx) = invalid;
            numCorrespondences = numCorrespondences + 1;
            indices(numCorrespondences) = astar;
        end
    end

      
        
        solEig = abs(V(:,1));
        %solEig = solEig(sol);
        [solEig, indx] = sort(solEig, 'descend');
        figure;
        subplot(2,1,1), plot(solEig);
        title('principal eigen vector');

        n = length(solEig) -1;
        difference = zeros(n-1,1);
        threshIndx = numCorrespondences;
        
        threshold = mean(diff(solEig));
        for i = 2:n
            difference(i-1) = solEig(i-1)-solEig(i);
            if(i > 2 && abs(difference(i-1) - difference(i-2)) > threshold)
                threshIndx = i-1;
                break;
            end
        end
        subplot(2,1,2), plot(difference);
        figure;
        plot(diff(solEig));
        title('consecutive diff of eigen vector values');
        
        l = l(indx,:);
        indices = indices(indx);
        % final correspondences based on solution indicator vector
        correspondences = l(1:threshIndx,:);  
        indices = indices(1:threshIndx);
        
end