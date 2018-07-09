function neighbourhood = find_neighbourhood(A, node, neighbourLevel)
%    neighbourLevel = 1;
%    A = [0 1 1 0 0 1 0 1 1 0 0 0
%         1 0 1 0 0 0 0 0 1 0 0 0
%         1 1 0 1 0 1 0 0 0 0 0 0
%         0 0 1 0 1 1 0 0 0 0 0 0
%         0 0 0 1 0 1 1 0 0 0 0 0
%         1 0 1 1 1 0 1 1 0 0 0 0
%         0 0 0 0 1 1 0 1 0 0 1 1
%         1 0 0 0 0 1 1 0 1 1 1 0
%         1 1 0 0 0 0 0 1 0 1 0 0
%         0 0 0 0 0 0 0 1 1 0 1 0
%         0 0 0 0 0 0 1 0 0 1 0 1
%         0 0 0 0 0 0 1 0 0 0 1 0];
 
    n = size(A,1);
    %node = randi([1 n],1,1)
    seenNodes = zeros(n,1);
    seenNodes(node) = 1;
    list = zeros(n,1);
    list(1) = node;
    level = zeros(n,1);
    last = 1;
    level(node) = 1;
    flag = 0;
    while ~all(seenNodes == 1)%length(find(seenNodes == 0)) > 0
        
        src = list(1);
        if (level(src) + 1) > neighbourLevel + 1
            'breaking';
            flag = 1;
            break;
        end
        if flag == 1
            break;
        end
        list(1) = [];
        last = last - 1;
        neighbours = find(A(src,:) > 0);
        %parent(neighbours) = src;
        for k = neighbours
            if seenNodes(k) ~= 1
                level(k) = level(src) + 1;
                last = last + 1;
                list(last) = k;     
                seenNodes(k) = 1;
            end
        end
    end
    
    neighbourhood = find(level > 0);
end
