close all;
clear all;
clc;

%% readImage
fileName = 'graph2.png';
fprintf('Reading file       : %s\n',fileName);
img = imread(fileName);
[m,n] = size(img);
vc = 0;
for i = 1:m
    for j = 1:n
        if img(i,j)==0
            vc = vc +1;
        end
    end
end
fprintf('Number of vertices : %d\n',vc);
vertices = zeros(vc,2);
index = 1;
for i = 1:m
    for j = 1:n
        if img(i,j)== 0
            vertices(index,:) = [i,j];
            index = index + 1;
        end
    end
end

%% create adjacency matrix of full graph
const = 10;
adjMatrix = zeros(vc);
for i = 1:vc
    for j = 1:vc
        adjMatrix(i,j) = exp(-norm(vertices(j,:)-vertices(i,:))/const);
        if i==j
            adjMatrix(i,j) = 0;
        end
    end
end
fprintf('Distance Matrix Calculated.\n');
return;
