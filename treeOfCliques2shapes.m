clc;
close all;
clear;

% load shapes
shapeName = {'david2.mat', 'david3.mat', 'david8.mat','david5.mat','david6.mat', 'david7.mat'};
numShapes = 2;
data = load(char(shapeName(1)));
numNodes = zeros(numShapes,1);
for i = 1:numShapes
    data(i) = load(char(shapeName(i)));
    numNodes(i) = length(data(i).surface.X);
end


A = zeros(numNodes(1), numNodes(1), numShapes);
L = zeros(numNodes(1), numNodes(1), numShapes);
V = zeros(numNodes(1), numNodes(1)-1, numShapes);
E = zeros(numNodes(1)-1, numNodes(1)-1, numShapes);
for i = 1:numShapes
    adjName = ['Adj_' char(shapeName(i))];
    t = load(char(adjName));
    A(:,:,i) = t.A;
    lapName = ['lap_' char(shapeName(i))];
    t1 = load(char(lapName));
    L(:,:,i) = t1.L;
    V(:,:,i) = t1.V;
    %[L(:,:,i),V(:,:,i), E(:,:,i)] = laplacian(t.A);
end

disp('using random initialization for clustering....');
k = 7;
d = k;
v = zeros(size(V,1), d, numShapes);
embeddingData = zeros(size(V, 1), size(V, 2), numShapes);
for i = 1:numShapes
        v(:,:,i) = V(:,1:d,i);
        embeddingData(:,:,i) = V(:,:,i);
end


%gmm1 = fitgmdist(v(:,:,1), k, 'Options',statset('MaxIter', 500), 'Replicates', 500);

sn1 = char(shapeName(1));
sn1 = char(sn1(1:find(sn1 == '.')-1));
sn1 = [sn1 '_gauss_' num2str(k) '_num2.mat'];
gmm1 = load(sn1);
gmm1 = gmm1.gmm;

%mu : rows:gaussian cols:mean dim 3rd dim: shape number
mu = zeros(size(gmm1.mu,1), size(gmm1.mu,2), numShapes);
sigma = zeros(size(gmm1.sigma,1), size(gmm1.sigma,2), size(gmm1.sigma,3), numShapes);
label = zeros(numNodes(1), numShapes);
mu(:,:,1) = gmm1.mu;
sigma(:,:,:,1) = gmm1.sigma;
label(:,1) = gmm1.label;%cluster(gmm1, v(:,:,1));
disp(['Mixture models created for shape ' num2str(1)]);
k = length(unique(label(:,1)));



map = [1, 1, 0
    0, 1, 1
    1, 0, 1
    1, 0, 0
    0, 1, 0
    0, 0, 1
    0.5, 0.5, 0.5
    0.2, 0.2, 0.2
    0.1, 0.8, 0.4
    1, 1, 1
    0.8 0.5 0.5
    0.3 0.7 0.9
    0.6 0.6 1
    0.3 0.2 0.7
    0.7 0.1 0.1
    ];
labels = {'clust1','Clust2','clust3','clust4','clust5','clust6','clust7','clust8','clust9','clust10','clust11','clust12','clust13','clust14','clust15'};


%load the gmm clustering matrices
for i = 2:numShapes
    sni = char(shapeName(i));
    sni = char(sni(1:find(sni == '.')-1));
    sni = [sni '_gauss_' num2str(k) '_num2.mat'];
    gmmi = load(sni);
    gmmi = gmmi.gmm;

    %gmmi = fitgmdist(v(:,:,i), k, 'Options',statset('MaxIter',500), 'Replicates', 500);
    mu(:,:,i) = gmmi.mu;
    sigma(:,:,:,i) = gmmi.sigma;
    label(:,i) = gmmi.label;%cluster(gmmi, v(:,:,i));
    disp(['Mixture models created for shape ' num2str(i)]);
    tempK = length(unique(label(:,i)));
    if tempK < k
        k = tempK;
    end
end



% % scale gaussian mean to lie between -1 and 1
% for i=1:numShapes
%     %mean of all gaussians of shape i
%     temp = mu(:,:,i);
%     %subtract the min of gaussian means from all gaussians for mean to lie
%     %from 0
%     zerocen = temp - repmat(min(temp), size(temp,1), 1);
%     df = max(temp) - min(temp)
%     normcen = zeros(size(temp));
%     for j = 1:size(temp,2)
%         normcen(:,j) = zerocen(:,j) ./ df(j);
%     end
%     scalecen = normcen*2 + repmat(-1, size(temp,1), size(temp,2));
%     mu(:,:,i) = scalecen;
% end
% mu

%display the gmm clustering
figure;
colormap(map(1:k,:));
trisurf(data(1).surface.TRIV, data(1).surface.X, data(1).surface.Y, data(1).surface.Z, label(:,1));

%trisurf(d1.surface.TRIV, d1.surface.X, d1.surface.Y, d1.surface.Z, newLabel1, 'Edgecolor', 'none', 'facecolor', 'interp' );
lcolorbar(labels(1:k));
translate = 150;
for i = 2:numShapes
    hold on;
    disp(['displaying surface ' num2str(i)]);
    trisurf(data(i).surface.TRIV, data(i).surface.X + repmat([translate], size(data(i).surface.X, 1), 1), data(i).surface.Y, data(i).surface.Z, label(:,i));
    translate = translate + 100;
end



%% association graph between shapes

% association between pairwise shapes
graphNodes = k*k;
newGraphNodes = graphNodes*(numShapes-1);
M = zeros(newGraphNodes, newGraphNodes);
l = zeros(graphNodes, 2);
j = 0;
dist = zeros(newGraphNodes, newGraphNodes);
degNodes = zeros(numShapes-1, graphNodes);
powerMij = zeros(graphNodes, graphNodes, numShapes-1);
normalizedMij = zeros(graphNodes, graphNodes, numShapes -1);
power = 5;
for i = 1: numShapes-1
    %[M, l]  = associationGraphSequence3(mu(:,:,i:i+1), sigma(:,:,:,i:i+1), 2, k, 1);
    [Mij, lT, distij] = associationGraphSequence3(mu(:,:,i:i+1), sigma(:,:,:,i:i+1), 2, k, 1);
    
    
    % degree of association graph;
    degNodes(i,:) = sum(Mij, 2);
    
    %power of association graph
    
    for p = 1: length(degNodes(i,:))
        normalizedMij(p,:,i) = Mij(p,:) ./ sum(Mij(p,:));
    end
    powerMij(:,:,i) = mpower(normalizedMij(:,:,i),power);
    
    figure;
    imagesc(Mij);
    %title('unnormalized Mij');
    
    %Mij = exp(-(distij.*distij)/std(distij(:)));
    %Mij(logical(eye(size(Mij)))) = 0;
    
    figure;
    imagesc(Mij);
    title('Normalized Mij');
    
    M(j+1:j+graphNodes, j+1:j+graphNodes) = Mij;
    l(j+1:j+graphNodes,:) = lT + (k*(i-1));
    dist(j+1:j+graphNodes, j+1:j+graphNodes) = distij;
    j = j + graphNodes;
end
figure;
imagesc(M);
title('Affinity matrix M');


%% Solution for the matching with Association Graph
%correspondences between shape1 and shape2
[correspondences] = correspondenceFromAssociation2(M, l, k*(numShapes-1));
cliqueWeight = findCliqueWeight(M, correspondences, k);
drawCorrespondences(data, label, labels, correspondences, k, numShapes, numNodes, cliqueWeight);

%% find the symmetric points
[symVertices, interestPoints, features] = findSymmetries(data, shapeName, numShapes);

%% find Association graph for symmetric points
% Calc weights of the graph using associationGraphSequence3
% (func from ShapeMatching)
associationSymm = containers.Map;
distSymm = containers.Map;
lSymm = containers.Map;
for i = 1:numShapes
    embeddingData = V(:,:,i);
    vertices = symVertices(num2str(i));
    [associationSymm(num2str(i)), distSymm(num2str(i)), lSymm(num2str(i))] = associationGraphSymmetryEuclidean(vertices, embeddingData, features(num2str(i)), interestPoints(num2str(i)));
    %figure;
    %imagesc(associationSymm(num2str(i)));
    %title(['Symmetry Association graph for shape ' num2str(i)]);
end
%}

%% find best symmetries from symmetric Association graph
symmetry = containers.Map;
symmIndices = containers.Map;
lIndx = containers.Map;
for i = 1:numShapes
    numVert = size(symVertices(num2str(i)),1);
    %[symmetry(num2str(i)), symmIndices(num2str(i))] = correspondenceFromAssociation2(associationSymm(num2str(i)), lSymm(num2str(i)), numVert);
    symmetry(num2str(i)) = correspondenceFromAssociationSymm(associationSymm(num2str(i)), lSymm(num2str(i)));
end


%find the gaussians of the symmetric points
%indx = keys(symVertices);
symIntraGauss = containers.Map;
symmPoints = containers.Map;

[sx,sy,sz] = sphere;
sphere_scale = 5;
sx = sphere_scale*sx;
sy = sphere_scale*sy;
sz = sphere_scale*sz;
for i = 1:numShapes
        % symmetric points
        %vertices = symVertices(num2str(i));
        vertices = symmetry(num2str(i));
        
        % take 60% of the vertices returned
%         symmthresh = 0.3;
%         ind = ceil(size(symmetry(num2str(i)), 1)* symmthresh);
%         vertices = vertices(1:ind,:);
%         
        % gaussian labels of all points in shape i 
        gausslabel = label(:,i);
        
        % variable to store the gaussian numbers of all vertex points in
        % shape i
        temp = zeros(length(vertices),2);
        tempSymmPoints = zeros(length(vertices),2);
        cntr = 1;
        % find gaussian labels of each row in vertices of shape i 
        for j = 1:size(vertices,1)
            % symmetry useful only between different gaussians
            if (label(vertices(j,1),i) ~= label(vertices(j,2),i))
                temp(cntr,:) = [label(vertices(j,1),i) label(vertices(j,2),i)];
                tempSymmPoints (cntr,:) = [vertices(j,1) vertices(j,2)];
                cntr = cntr + 1;
                
                color = label(:,i);
                color(vertices(j, 1)) = max(gausslabel) + 1;
                color(vertices(j, 2)) = max(gausslabel) + 1;
                %{
                tx1 = data(i).surface.X(vertices(j,1));
                ty1 = data(i).surface.Y(vertices(j,1));
                tz1 = data(i).surface.Z(vertices(j,1));
                tx2 = data(i).surface.X(vertices(j,2));
                ty2 = data(i).surface.Y(vertices(j,2));
                tz2 = data(i).surface.Z(vertices(j,2));
                %}
                figure;
                trisurf(data(i).surface.TRIV, data(i).surface.X, data(i).surface.Y, data(i).surface.Z, color);
                %{
                hold on;
               surf(sx + tx1, sy + ty1, sz + tz1);
                hold on;
                surf(sx + tx2, sy + ty2, sz + tz2);
                %}
                title([num2str(label(vertices(j,1),i)) ' ' num2str(label(vertices(j,2),i))]);
                shading interp;
                
            end
        end
        
        % take only unique correspondences
        [temp, ia] = unique(temp, 'rows');
        tempSymmPoints  = tempSymmPoints (ia,:);
        
        %remove any invalid(zero) correspondences
        in = find(temp(:,1) ~= 0);
        temp = temp(in,:);
        tempSymmPoints  = tempSymmPoints (in,:);
        in = find(temp(:,2) ~= 0);
        temp = temp(in,:);
        tempSymmPoints  = tempSymmPoints (in,:);
        
        % store all gaussian numbers in the map
        symIntraGauss(num2str(i)) = temp;    
        
        %store all symmetric points in the map
        symmPoints(num2str(i)) = tempSymmPoints;
end



% store all solutions
global symmSols;
symmSols = containers.Map;
cntr = num2str(1);
symmSols(cntr) = correspondences;

% find new solution for each shape
for i = 1:numShapes
    symmGauss = symIntraGauss(num2str(i));
    findSymmSol(symmGauss, i); 
end


%find unique sols from the solution set
keys = symmSols.keys;
uniqueSols = zeros(size(correspondences,1), size(correspondences, 2), length(keys));
uniqueSols(:,:,1) = symmSols(char(keys(1)));
uniqueNum = 1;
for i = 2:length(keys)
    insert = true;
    for j = 1:uniqueNum
        if(isequal(sortrows(symmSols(char(keys(i)))), sortrows(uniqueSols(:,:,j))))
           insert = false;
        end
    end
    if(insert)
        uniqueNum = uniqueNum + 1;
        uniqueSols(:,:, uniqueNum) = symmSols(char(keys(i)));
    end
end
uniqueSols = uniqueSols(:,:, 1:uniqueNum);

%check for correct sol
correctSol01 = sortrows([ 1 1; 2 7; 3 5; 4 4; 5 2; 6 3; 7 6]);
correctSol12 = sortrows([1 2; 2 3; 3 6; 4 1; 5 5; 6 7; 7 4]);
for i = 1:uniqueNum
    
    sol = uniqueSols(:,:,i);
    cliqueWeight = findCliqueWeight(M, uniqueSols(:,:,i), k);
    color1 = zeros(numNodes(1),1);
    color2 = zeros(numNodes(1),1);
    for j = 1:k
        label1 = find(label(:,1) == sol(j,1));
        label2 = find(label(:,2) == sol(j,2));
        color1(label1) = j;
        color2(label2) = j;
       
    end
    
    sol = sortrows(uniqueSols(:,:,i));
    if(sum(sum(sol == correctSol01)) == numel(correctSol01))
        disp(['correct solution found clique weight ' num2str(cliqueWeight)]);
    else
        disp(['clique weight ' num2str(cliqueWeight)]);
    end
    
     figure;
     trisurf(data(1).surface.TRIV, data(1).surface.X, data(1).surface.Y, data(1).surface.Z, color1);
     hold on;
     trisurf(data(2).surface.TRIV, data(2).surface.X  + repmat([150], numNodes(1),1), data(2).surface.Y, data(2).surface.Z, color2);
     
     title(['Weight of Clique ' num2str(cliqueWeight)]);
     
end
%} 