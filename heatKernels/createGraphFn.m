function [adjMatrix] = createGraphFn(fileName,k,type,vc)

%% readImage

img = imread(fileName);
[m,n] = size(img);

%number of vertices
vc = length(find(img == 0));

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
adjMatrix = zeros(vc,vc);

%{
idx = find(img == 0);
[X,Y] = ind2sub([vc 2], idx);
vert = zeros(vc,2);
vert(idx,:) = [X,Y]; 
isequal([X,Y], vertices)
vert'
'**************'
vertices'
%}


%% create adjacency matrix of full graph
adjMatrix = zeros(vc,vc);
if type == 1
   sig = 0.5;
   pd = pdist2(vertices,vertices,'euclidean');
   pd = pd/(max(pd(:))-min(pd(:)));
   adjMatrix1 = exp(-pd.^2/(2*sig^2));
   adjMatrix1(logical(eye(size(adjMatrix1)))) = 0;
   adjMatrix = adjMatrix1;
elseif type == 2
   sig = 0.5;
   pd = pdist2(vertices,vertices,'cityblock');
   pd = pd/(max(pd(:))-min(pd(:)));
   adjMatrix1 = exp(-pd/(2*sig^2));
   adjMatrix1(logical(eye(size(adjMatrix1)))) = 0;
   adjMatrix = adjMatrix1;
end
%%euclidean dist



fprintf('Distance Matrix Calculated.\n');
if k==0
    return;
else
    
    % considering top k neighbours.
    adjMat2 = zeros(size(adjMatrix));
    for i = 1:vc
        j = 0;
        while j < k
            dat = find(adjMatrix(i,:) == max(adjMatrix(i,:)));
            if length(dat) > k
                dat = dat(1:k);
                j = j + k;
            else
               j = j + length(dat);
            end
               adjMat2(i,dat) = adjMatrix(i,dat);
               adjMatrix(i,dat) = 0;
        end
    end
    fprintf('Adjacent Matrix:\n');
    adjMatrix=adjMat2;
    
end

plotGraphEdges(vertices, adjMatrix);

end