function edges =  find_edges(surface)

    triangles = surface.TRIV;
    
    % number of nodes
    n = length(surface.X);
    
    edges = zeros(n*n,2);

    unqeEdges = unique(sort(triangles,2), 'rows');
        
    temp = unique(sort(unqeEdges(:,1:2),2), 'rows');
    edges(1:length(temp),:) = temp;
    l = length(temp);
    
    temp = unique(sort(unqeEdges(:,2:3),2), 'rows');
    edges(l+1:l+length(temp),:) = temp;
    l = length(temp) + l;
    
    temp = unique(sort(unqeEdges(:,3:1),2), 'rows');
    edges(l+1:l+length(temp),:) = temp;
    l = length(temp) + l;
    
    edges = unique(edges(1:l,:),'rows');
    
%     temp = unique(sort(triangles(:,1:2),2), 'rows');
%     edges(1:length(temp),:) = temp;
%     l = length(temp);
%     
%     temp = unique(sort(triangles(:,2:3),2), 'rows');
%     edges(l+1:l+length(temp),:) = temp;
%     l = length(temp) + l;
%     
%     temp = unique(sort(triangles(:,3:1),2), 'rows');
%     edges(l+1:l+length(temp),:) = temp;
%     l = length(temp) + l;
%     
%     edges = unique(edges(1:l,:),'rows');
%     
end