%% Computes Adjacency Matrix from traingluar meshes

function A = adjacency_matrix(surface)
    
    % number of nodes
    n = length(surface.X);
    
    edges = find_edges(surface);
  
    A = zeros(n,n);
    u = [surface.X(edges(:,1)) surface.Y(edges(:,1)) surface.Z(edges(:,1))];
    v = [surface.X(edges(:,2)) surface.Y(edges(:,2)) surface.Z(edges(:,2))];
    sigma = mean(sqrt(sum((u-v).^2,2)));
    edges_weight =  exp(-sqrt(sum((u-v).^2,2))/sigma);
    
    for i = 1:length(edges)
        A(edges(i,1), edges(i,2)) = edges_weight(i);
        A(edges(i,2), edges(i,1)) = edges_weight(i);
    end
    
end
