function [  ] = plotGraph( idx, fileName )
%PLOTGRAPH Plot graph as per the fileName
%   Give color corresponding to the idx.

%% readImage
fprintf('plotGraph :: Reading file       : %s\n',fileName);
img = imread(fileName);
[m,n] = size(img);
vc = length(idx);
fprintf('plotGraph :: Number of vertices : %d\n',vc);

vertices = zeros(vc,2);
index = 1;

for i = 1:m
    for j = 1:n
        if img(i,j)== 0
            vertices(index,:) = [(m-i),j];
            index = index + 1;
        end
    end
end

%% form the color matrix
colVertex=zeros(vc,3);
for i = 1:vc
    if idx(i) == 1
        colVertex(i,:)=[0,1,1];
    elseif idx(i) == 2
        colVertex(i,:)=[1,1,0];
    elseif idx(i) == 3
        colVertex(i,:)=[1,0,1];
    else
        colVertex(i,:)=[1,0,0];
    end
end
figure('name','clusters formed');
scatter(vertices(:,2),vertices(:,1),10,colVertex,'filled');

axis([0 m 0 n]);

end

