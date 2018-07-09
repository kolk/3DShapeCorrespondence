function [symm_graph_vertices, interestPoints, probDistShapes] = findSymmetries(shapes, shapeName, numShapes)
    %{
Program to create a graph that records the symmetry of the object using
Heat Kernels
%}

% clc;
% close all;
% clear all;
% shapeName = {'david5.mat', 'david1.mat', 'david8.mat','david5.mat','david6.mat', 'david7.mat'};
% numShapes = 1;
% shapes = load(char(shapeName(1)));
% for i = 1:numShapes
%     shapes(i) = load(char(shapeName(i)));
%     
% end

%% Find Symmetric points
numNodes = length(shapes(1).surface.X);
A = zeros(numNodes, numNodes, numShapes);
L = zeros(numNodes, numNodes, numShapes);
D = zeros(numNodes, numNodes, numShapes);
H = zeros(numNodes, numNodes, numShapes);
eig_vec = zeros(numNodes, numNodes-1, numShapes);
eig_val = zeros(numNodes-1,numNodes-1, numShapes);
symm_graph_vertices = containers.Map;
interestPoints = containers.Map;
probDistShapes = containers.Map;


for i = 1:numShapes
    adjName = ['Adj_' char(shapeName(i))];
    t = load(char(adjName));
    A(:,:,i) = t.A;
    lapName = ['lap_' char(shapeName(i))];
    t1 = load(char(lapName));
    L(:,:,i) = t1.L;
    eig_vec(:,:,i) = t1.V;
    eig_val(:,:,i) = t1.E;

    ev = diag(t1.E);
    lambda = ev(numNodes(1)*0.01);
    heatDiffVal = 0.05; 
    scale = -log(heatDiffVal)/lambda
    %scale = 50;
  
    % Calulating Heat Kernel Matrix
    D = diag(exp(-diag(eig_val(:,:,i)) * scale));
    H(:,:,i) = eig_vec(:,:,i) * D * eig_vec(:,:,i)';

    % Auto-diffusion function
    autoDiff = diag(H(:,:,i));

    % find extrema of ADF
    [sortedAdf, nodes] = sort(autoDiff, 'descend');

    neighbourhood = 4;
    extm = getExtrema(A(:,:,i), nodes, neighbourhood);
    ip = find(extm == 1);
    
    % store interst points for shape i
    interestPoints(num2str(i)) = ip;
    
    epsilone=10^-5;
    
    % The rows of the heat kernel corresponding to the interest points
    H_ip = H(ip, :,i);


    % label the interest points 

    color = zeros(numNodes,1);
    map = [0 0 1;
       1 0 1;
       1 1 0
           ];
    color(ip) = 1;

    [sx,sy,sz] = sphere;
    sphere_scale = 2;
    sx = sphere_scale*sx;
    sy = sphere_scale*sy;
    sz = sphere_scale*sz;
    
    figure;
    colormap(map);
    trisurf(shapes(i).surface.TRIV, shapes(i).surface.X, shapes(i).surface.Y, shapes(i).surface.Z, color);
    for j = 1:length(ip)
        tx1 = shapes(i).surface.X(ip(j));
        ty1 = shapes(i).surface.Y(ip(j));
        tz1 = shapes(i).surface.Z(ip(j));
    hold on;
    surf(sx + tx1, sy + ty1, sz + tz1);
    end
    
    shading interp;
    camlight;
    title('Interest points ');

    % compute prob density metric
    y = H_ip(1,H_ip(1,:)>epsilone);
    %y = H_ip(1,:);
    [pd xi] = ksdensity(y);
    probDist = zeros(size(H_ip,1), length(pd)); 
    probDistXi = zeros(size(H_ip,1), length(xi));
       
   for j = 1:size(ip, 1)
        y = H_ip(j,H_ip(j,:)>epsilone);
        %y = H_ip(j,:);
        [pd, xi] = ksdensity(y);
        probDist(j,:) = (pd - mean(pd))/std(pd);
       
        probDistXi(j,:) = xi;
    end

    % compute the pairwise heatkernel prob densities with the bhattacharya
    % dist metric
    comp = zeros(size(ip,1), size(ip, 1)); 
    for j = 1:size(ip, 1)
        for k = j:size(ip, 1)
            comp(j,k) = sqrt(sum( probDist(j,:) .* probDist(k,:)));
%             H_ip(j,:) = (H_ip(j,:) - mean(H_ip(j,:)))/std(H_ip(j,:));
%             H_ip(k,:) = (H_ip(k,:) - mean(H_ip(k,:)))/std(H_ip(k,:));
%             comp(j,k) = sqrt(sum(H_ip(j,:) - H_ip(k,:).^2));
            %comp(j,k) = sqrt(sum( probDist(j,:) - probDist(k,:)));
            comp(k,j) = comp(j,k);
        end
    end

    % store the prob densities
    probDistShapes(num2str(i)) = probDist;
    
    % Each pair of interst points form a vertex if their bhat dist is within epsilon difference
    vertices = [];
    dist = [];
    for j = 1:size(ip, 1)
        %indices = find(comp(cntr,:) < epsilon);
        [d, temp] = sort(comp(j, :), 'descend'); 
   
        %eliminate comparison with the same node
        temp = temp(find(temp ~= j));

        for k = 1:2
            vertices = [vertices; j temp(k)];
            dist = [dist; d(k)];
        end
    end
    
    [~, indx] = sort(dist);
    vertices = vertices(indx,:);
    
    %symm_rels = vertices;
    symm_rels = [];
    s = 1;
    % Finding symmetric relations
    for j = 1:size(vertices, 1);
        idx = vertices(j, 1);
        val = vertices(j, 2);   
    
        % add node to graph only if it has a symmetric counterpart
        if sum(ismember(vertices(find(vertices(:, 1) == val), 2), idx)) ~= 0
            % if node is already present as part of a symmetric relation, dont
            % add
            if(length(find(ismember(symm_rels, [val idx])) > 0) == 0)
                symm_rels(s, :) = [idx val];
                s = s + 1;
            end
        end
    end

    % Corresponding graph vertices
    graph_ver = zeros(size(symm_rels,1), 2);
    graph_ver(:, 1) = ip(symm_rels(:, 1));  
    graph_ver(:, 2) = ip(symm_rels(:, 2));
    
    symm_graph_vertices(num2str(i)) = graph_ver;
    
    
    %color points of a vertex with same color
   %{
    for j = 1:size(symm_rels,1)
        col1 = zeros(numNodes,1);
        col1(graph_ver(j, 1)) = i;
        col1(graph_ver(j, 2)) = i;
        figure; 
        trisurf(shapes(i).surface.TRIV, shapes(i).surface.X, shapes(i).surface.Y, shapes(i).surface.Z, col1);
    
        %index of vertex in comp
        indx1 = find(ip == ip(symm_rels(j,1)))
        indx2 = find(ip == ip(symm_rels(j,2)))  
        
        title(['comp ' num2str(comp(indx1, indx2))]);
              
        
        tx1 = shapes(i).surface.X(graph_ver(j,1));
        ty1 = shapes(i).surface.Y(graph_ver(j,1));
        tz1 = shapes(i).surface.Z(graph_ver(j,1));
        tx2 = shapes(i).surface.X(graph_ver(j,2));
        ty2 = shapes(i).surface.Y(graph_ver(j,2));
        tz2 = shapes(i).surface.Z(graph_ver(j,2));
        hold on;
        surf(sx + tx1, sy + ty1, sz + tz1);
        hold on;
        surf(sx + tx2, sy + ty2, sz + tz2);
        shading interp;
    end
    %}
end

end

