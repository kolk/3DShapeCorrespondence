function [L,V,E] = laplacian(A)
    %% construct degree Matrix
    deg = sum(A,2);
    D = diag(deg);
    
    %% construct laplacian and eigenVec
    
    
    %unnormalized laplacian
    L = D - A;
    
    %normalized laplacian
    %deg(deg == 0) = 0.000001;
    %D = diag(1./(deg.^0.5));
    
    %L = inv(D).^0.5 * L * inv(D).^0.5;
    
    [V,E] = eig(L); % V - Eigenvectors E - Eigenvalues
    [E, idx] = sort(diag(E));
    V = V(:,idx);
    V=V(:,2:end);
    E = diag(E(2:end));

end