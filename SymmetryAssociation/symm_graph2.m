%{
Program to create a graph that records the symmetry of the object using
Heat Kernels
%}
%function symm_graph2(shapes, numNodes, A)
clc;
close all;
clear all;

shapeName = {'david0.mat', 'david1.mat', 'david8.mat','david5.mat','david6.mat', 'david7.mat'};
numShapes = 2;
shapes = load(char(shapeName(1)));
numNodes = zeros(numShapes,1);
for i = 1:numShapes
    shapes(i) = load(char(shapeName(i)));
    numNodes(i) = length(shapes(i).surface.X);
end

% same size shapes
numNodes = numNodes(i);

A = zeros(numNodes, numNodes, numShapes);

L = zeros(numNodes, numNodes, numShapes);
D = zeros(numNodes, numNodes, numShapes);
H = zeros(numNodes, numNodes, numShapes);
eig_vec = zeros(numNodes, numNodes-1, numShapes);
eig_val = zeros(numNodes-1,numNodes-1, numShapes);
symm_graph_vertices = containers.Map;
scale = 20;%7.2316;
for i = 1:numShapes
    %adjacency ma
    adjName = ['Adj_' char(shapeName(i))];
    t = load(char(adjName));
    A(:,:,i) = t.A;
    [L(:,:,i),eig_vec(:,:,i),eig_val(:,:,i)] = laplacian(t.A);
    
    % Calulating Heat Kernel Matrix
    D = diag(exp(-diag(eig_val(:,:,i)) * scale));
    H(:,:,i) = eig_vec(:,:,i) * D * eig_vec(:,:,i)';

    % Auto-diffusion function
    autoDiff = diag(H(:,:,i));

    % find extrema of ADF
    [sortedAdf, nodes] = sort(autoDiff, 'descend');

    neighbourhood = 2;
    extm = getExtrema(A(:,:,i), nodes, neighbourhood);
    ip = find(extm == 1);
    
    epsilone=10^-5;
    % The rows of the heat kernel corresponding to the interest points
    H_ip = H(ip, :,i);

    %data in embedding space
    dimension = 7;
    embeddingData = eig_vec(:,1:dimension,:);

    % label the interest points 

    color = zeros(numNodes,1);
    map = [0 0 1;
       1 0 1;
       1 1 0
           ];
    color(ip) = 1;

    figure;
    colormap(map);
    trisurf(shapes(i).surface.TRIV, shapes(i).surface.X, shapes(i).surface.Y, shapes(i).surface.Z, color);
    shading interp;
    camlight;
    title('Interest points ');

    y = H_ip(1,H_ip(1,:)>epsilone);
    [pd xi] = ksdensity(y);
    probDist = zeros(size(H_ip,1), length(pd)); 
    probDistXi = zeros(size(H_ip,1), length(xi));
    probDist(1,:) = pd;
    probDistXi(1,:) = xi;
    for j = 2:size(ip, 1)
        y = H_ip(j,H_ip(j,:)>epsilone);
        [pd, xi] = ksdensity(y);
        probDist(j,:) = pd;
        probDistXi(j,:) = xi;
    end


    comp = zeros(size(ip,1), size(ip, 1)); 
    for j = 1:size(ip, 1)
        for k = j:size(ip, 1)
            pd1 = (probDist(j,:) - mean(probDist(j,:)))/std(probDist(j,:));
            pd2 = (probDist(k,:) - mean(probDist(k,:)))/std(probDist(k,:));
            comp(j,k) = sqrt(sum( pd1 .* pd2));
            comp(k,j) = comp(j,k);
        end
    end


    % Each pair of interst points form a vertex if their bhat dist is within epsilon difference
    vertices = [];
    epsilon = 0.01;
    for j = 1:size(ip, 1)
        %indices = find(comp(cntr,:) < epsilon);
        [~, temp] = sort(comp(j, :), 'descend'); 
   
        %eliminate comparison with the same node
        temp = temp(find(temp ~= j));

        vertices = [vertices; j temp(1)];
    end
    
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
        shading interp;
    
    end
    %}
end
%}

%% find Association graph
% Calc weights of the graph using associationGraphSequence3
% (func from ShapeMatching)

M = containers.Map;
dist = containers.Map;
l = containers.Map;
for i = 1:numShapes
    [M(num2str(i)), dist(num2str(i)), l(num2str(i))] = associationGraphSymmetry(symm_graph_vertices(num2str(i)), embeddingData(:,:,i));
end
%% find correspondence from Association graph
correspondences = containers.Map;
indices = containers.Map;
lIndx = containers.Map;
for i = 1:numShapes
    correspondences(num2str(i)) = correspondenceFromAssociation2(M(num2str(i)), l(num2str(i)), n);
end


%end