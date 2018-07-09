clear;
 close all;

% load shape1 and shape2
 d1 = load('david4.mat');

 
% number of nodes
n1 = length(d1.surface.X);

    
%adjacency matrix
A1 = adjacency_matrix(d1.surface);


% laplacian and eigen decomposition
[L1,V1,E1] = laplacian(A1); %V eigen vectors E- eigen values


[~,n] = size(V1);
v = diag(E1);
spec = {'r-'; 'b-'; 'g-'; 'y-'; 'c-'; 'k-'; 'm-'; 'r--'; 'b--'; 'g--'};

k = 1;
figure;
t = 2;
%for t = 1:50:300
    xlabel('lambda');
    ylabel('exp(-lamda*t)');
    y = zeros(1,n);
    for i = 1:n
        y(i) = exp(v(i)*t);
    end
    
    spec(k)
    y
    plot(1:n, y, char(spec(k)));
    k = k+1;
    
    hold on;
%end    
    % construct heat kernel
% H1 = V1*diag(diag(exp(-E1*t)))*V1';

%randomly choose a heat scource
% source1 = randi([1 n1],1,1);
% 
% % intial heat distribution
% hiDistr1 = zeros(n1,1);
% hiDistr1(source1) = 1;
% 
% % final heat distribution
% hiDistr1 = H1*hiDistr1;
% hiDistr1 = hiDistr1/(max(hiDistr1)-min(hiDistr1));

% figure;
% colormap('jet');
% trisurf(d1.surface.TRIV, d1.surface.X, d1.surface.Y, d1.surface.Z, hiDistr1);
% shading interp;
% camlight;
% title(['scale parameter  ' num2str(t)]);
% shading interp;
% caxis([min(hiDistr1) max(hiDistr1)]);
% colorbar



%     % Auto-diffusion function
%     autoDiff1 = diag(H1);
% 
%     % find extrema of ADF
%     [sortedAdf1, nodes1] = sort(autoDiff1, 'descend');
% 
%     extm1 = getExtrema(A1, nodes1, 10);
% 
%     ind = find(extm1 == 1);
%     extremaClr1 = zeros(n1,1);
%     extremaClr1(ind) = 1;
%     size(ind)
%     figure;
%     colormap('jet');
%     trisurf(d1.surface.TRIV, d1.surface.X, d1.surface.Y, d1.surface.Z, extremaClr1);
%     shading interp;
%     %camlight;
%     title(['Extrema for scale parameter  ' num2str(t)]);

%legend('t 1', 't 3', 't ', 't1000');


