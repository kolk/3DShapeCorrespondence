function extm = getExtrema(A, nodes, neighbourhood)
    j = 1;
    n = size(A,1);
    seen = zeros(n,1);
    extm = zeros(n,1);
    while(~all(seen == 1))
        % if local neighbourhood of node is seen, move to next node
        node = nodes(j);
        if seen(node) ~= 1
            ngh = find_neighbourhood(A, node, neighbourhood);
            if all(seen(ngh) == 0)
                extm(node) = 1;   
                
            end
            seen(node) = 1;
            seen(ngh) = 1;
            %if(~ismember(ngh, seen))

        end
        j = j + 1;
    end

end