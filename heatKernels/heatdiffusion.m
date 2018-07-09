clear all;
close all;

%% create graph
connectivity=15;
fileName= 'graph4.png';
type = 1;

img = imread(fileName);
[m,n] = size(img);
vc = length(find(img == 0));
fprintf('plotGraph :: Number of vertices : %d\n',vc);

vertices = zeros(vc,2);
index = 1;

for i = 1:m
    for j = 1:n
        if img(i,j)== 0
            vertices(index,:) = [i,j];%[(m-i),j];
            index = index + 1;
        end
    end
end

adjMatrix = createGraphFn(fileName, connectivity, type,vc);


%% construct degree Matrix
degMatrix = diag(sum(adjMatrix,2));

%% construct laplacian
lapMatrix = degMatrix - adjMatrix;
 
%% construct hear kernel
t = 200;
scaledLap = -lapMatrix .* t;
H = expm(scaledLap);

%% iterative heat diffusion

%intial heat distribution
hiDistr = zeros(vc,1);

%randomly choose a heat scource
source = randi([1 vc],1,1); 

%set the distribution of source to 1
hiDistr(source) = 1;

%iteratively perform heat diffusion until convergence
prev = 0;
 while abs(prev - hiDistr(source)) > 0.00001
       prev = hiDistr(source);
       hiDistr = H*hiDistr;     
 end
    
 %% find clusters

k = 2;
[idx,C] = kmeans(hiDistr,k);
clust1 = find(idx == 1);
clust2 = find(idx == 2);

%plot the data
colVertex=zeros(vc,3);
colVertex(clust1,:) = repmat([1,0,0],[size(clust1,1) 1]);
colVertex(clust2,:) = repmat([0,0,1],[size(clust2,1) 1]);
figure;
scatter(vertices(:,2),vertices(:,1),10,colVertex,'filled');   