first_run=1;
if first_run==1
    clear;
    close all;
    
    %load cat0.mat
    load michael0.mat
     
    % number of nodes
    n = length(surface.X);
    
    % Adjacency Matrix
    A = adjacency_matrix(surface);
     
    
    % construct laplacian and eigenVec
    [L,V,E] = laplacian(A);
end

close all;

%randomly choose a heat scource
source = randi([1 n],1,1);
    
scale = [100:2000:10000];
for t = 1:length(scale)
    %% construct heat kernel
    H = V*diag(diag(exp(-E*scale(t))))*V';
    
    %% heat distribution
    
    %intial heat distribution
    hiDistr = zeros(n,1);
    hiDistr(source) = 1;
  
    %final heat distribution
    hiDistr = H*hiDistr;
    hiDistr=hiDistr/(max(hiDistr)-min(hiDistr));
    
    figure;
    colormap('jet');
    trisurf(surface.TRIV, surface.X, surface.Y, surface.Z, hiDistr);
    shading interp;
    camlight;
    title(['scale parameter  ' num2str(scale(t))]);
    shading interp;
    caxis([min(hiDistr) max(hiDistr)]);
    colorbar;
    %caxis([0 1])
     
 
    %% Auto-diffusion function
    autoDiff = diag(H);
    
    %% find extrema of ADF
    
    [sortedAdf,nodes] = sort(autoDiff, 'descend');
    %[h, node] = max(autoDiff);
    %node = 2892;
    seen = zeros(n,1);
    extm = zeros(n,1);
    
j = 1;
while(~all(seen == 1))
    % if local neighbourhood of node is seen, move to next node
    node = nodes(j);
    if seen(node) ~= 1
        ngh = find_neighbourhood(A, node, 2);
        if all(seen(ngh) == 0)
            extm(node) = 1;   
            disp('extrema assignment');
        end
        seen(node) = 1;
        seen(ngh) = 1;
        %if(~ismember(ngh, seen))
       
    end
    j = j + 1;
 end

    extremaClr = zeros(n,1);
    extremaClr(find(extm==1)) = 1;
   
    figure;
    colormap('jet');
    trisurf(surface.TRIV, surface.X, surface.Y, surface.Z, extremaClr);
    shading interp;
    %camlight;
    title(['Extrema for scale parameter  ' num2str(scale(t))]);
    
   
    
end
