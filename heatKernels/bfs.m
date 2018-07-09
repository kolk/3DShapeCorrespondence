function neighbourhood = bfs(A)
    n = size(A,1);
    node = randi([1 n],1,1);
    seenNodes(node) = 1;
    list = zeros(n,1);
    list(1) = node;
    level = zeros(n,1);
    last = 1;
    level(node) = 1;
    while length(find(seenNodes == 0)) > 0
        
        src = list(1);
        if (level(src) + 1) == neighbourLevel
            break;
        end
        
        list(1) = [];
        neighbours = find(A(src,:) > 0);
        for k = neighbours
            if seenNodes(k) ~= 1
                level(k) = level(src) + 1;
                last = last + 1;
                list(last) = k;               
            end
        end
    end
    neighbourhood = find(level > 0);
end