function [correspondences, indices] = correspondenceFromAssociation2(M, l, k)
    [V, D] = eig(M);
    D = diag(D);
    [D, ind] = sort(D,'descend'); 
    V = V(:,ind); 
    
    n = size(l,1);
    % solution vector
    sol = zeros(n,1);
    
    xstar = V(:,1)
    sort(xstar)
    
    pos = find(xstar > 0);
    neg = find(xstar < 0);
    doubleSol = false;
    
    
    if (length(pos) > 0 && length(neg) > 0)
        if(length(pos) > k && length(neg) > k)
            doubleSol = true;
            zeroCrossing = 0;
            sign = abs(xstar(1))/xstar(1);
            for i = 2:length(xstar)
                if sign ~= abs(xstar(i))/xstar(i)
                    zeroCrossing = i;
                    break;
                end
            end
        elseif(length(pos) > k)
            xstar = xstar(pos);
            l = l(pos,:);
        else
            xstar = xstar(neg);
            l = l(neg,:);
        end
    else
        xstar = abs(xstar);
    end
doubleSol
    
    [~, indx] = sort(xstar, 'descend');
    disp([l(indx,:) xstar(indx)])
    
    indices = [];
    if ~doubleSol
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
                indices = [indices; astar];
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
            end
        end

        % final correspondences based on solution indicator vector
        correspondences = l(sol>0,:)
       
    else
        return 
    end
    
end