%{
Program to create a graph that records the symmetry of the object using
Heat Kernels
%}

clc;
close all;
%clear all;

path = 'david2.mat';
shape = load(path);
n = size(shape.surface.X, 1);


%% Find Symmetric points
% Fartherst point interest points
%[pt, ip] = farthest_points([shape.surface.X shape.surface.Y shape.surface.Z], 200);
%indLeft = find((shape.surface.X > -80) & (shape.surface.X < -77.1) & (shape.surface.Y > -12.6) & (shape.surface.Y < 0) & (shape.surface.Z > 126) & (shape.surface.Z < 128.4)== 1);
%indLeft = indLeft(2:length(indLeft));
%indRight = find((shape.surface.X > 77) & (shape.surface.X < 80.3) & (shape.surface.Y >-12.6) & (shape.surface.Y < 0) & (shape.surface.Z > 126) & (shape.surface.Z < 128.4)== 1);
%indRight = indRight(1:2:length(indRight));
%ip = [indLeft; indRight; ip];



A = adjacency_matrix(shape.surface);
[Lap, eig_vec, eig_val] = laplacian(A);   

% Calulating Heat Kernel Matrix
t = 20;%7.2316;
D = diag(exp(-diag(eig_val) * t));
H = eig_vec * D * eig_vec';


% Auto-diffusion function
autoDiff = diag(H);

% find extrema of ADF
[sortedAdf, nodes] = sort(autoDiff, 'descend');

neighbourhood = 2;
extm = getExtrema(A, nodes, neighbourhood);
ip = find(extm==1);
%}
epsilone=10^-5;
% The rows of the heat kernel corresponding to the interest points
H_ip = H(ip, :);


% label the interest points 

color = zeros(n,1);
map = [0 0 1;
       1 0 1;
       1 1 0
       ];
color(ip) = 1;

figure;
colormap(map);
trisurf(shape.surface.TRIV, shape.surface.X, shape.surface.Y, shape.surface.Z, color);
shading interp;
camlight;
title('Interest points ');

%{
figure,
for i=1:length(ip)
    temp = H_ip(i,H_ip(i,:)>epesilone);
    %temp = smooth(H_ip(i,:));
    [F XI] = ksdensity(temp);
    subplot(8,5,i), plot(XI,F);
    ylim([0 1000]);
    title(num2str(i));
end
%}

y = H_ip(1,H_ip(1,:)>epsilone);
[pd xi] = ksdensity(y);
probDist = zeros(size(H_ip,1), length(pd));
probDistXi = zeros(size(H_ip,1), length(xi));
probDist(1,:) = pd;
probDistXi(1,:) = xi;
for i = 2:size(ip, 1)
    y = H_ip(i,H_ip(i,:)>epsilone);
    %y = H_ip(i,:);
    [pd, xi] = ksdensity(y);
    probDist(i,:) = pd;
    probDistXi(i,:) = xi;
end


comp = zeros(size(ip,1), size(ip, 1)); 
for i = 1:size(ip, 1)
    for j = i:size(ip, 1)
       pd1 = (probDist(i,:) - mean(probDist(i,:)))/std(probDist(i,:));
       pd2 = (probDist(j,:) - mean(probDist(j,:)))/std(probDist(j,:));
       comp(i,j) = sqrt(sum( pd1 .* pd2));
       comp(j,i) = comp(i,j);
    end
end


% Each pair of interst points form a vertex if their bhat dist is within epsilon difference
vertices = [];
epsilon = 0.01;
for i = 1:size(ip, 1)
   %indices = find(comp(cntr,:) < epsilon);
   [~, temp] = sort(comp(i, :), 'descend'); 
   
   %eliminate comparison with the same node
   temp = temp(find(temp ~= i));
   
   numNodes = 1;
   for indx = 1:numNodes
       vertices = [vertices; i temp(indx)];
   end    
end

%symm_rels = vertices;
symm_rels = [];
s = 1;
% Finding symmetric relations
for i = 1:size(vertices, 1);
    idx = vertices(i, 1);
    val = vertices(i, 2);   
    
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

%color points of a vertex with same color
for i = 1:size(symm_rels,1)
    col1 = zeros(n,1);
    col1(graph_ver(i, 1)) = i;
    col1(graph_ver(i, 2)) = i;
    figure; 
    colormap(map);
    trisurf(shape.surface.TRIV, shape.surface.X, shape.surface.Y, shape.surface.Z, col1);
    
    %index of vertex in comp
    indx1 = find(ip == ip(symm_rels(i,1)))
    indx2 = find(ip == ip(symm_rels(i,2)))  
    
    title(['comp ' num2str(comp(indx1, indx2))]);
    shading interp;
    camlight;
end


%% find Association graph
% Calc weights of the graph using associationGraphSequence3
% (func from ShapeMatching)
%[M, l] = associationGraphSymmetry(graph_ver, shape);


%% find correspondence from Association graph
%[correspondences, indices] = correspondenceFromAssociation2(M, l, k);
