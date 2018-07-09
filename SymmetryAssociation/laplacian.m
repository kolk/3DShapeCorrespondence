%% Calculating eigen values & vectors as well as Laplacian

function [L, eig_vec, eig_val] = laplacian(A)

    % construct degree Matrix
    D = diag(sum(A,2));
    
    % construct unnormalized laplacian and eigenVec
    L = D - A;
    [eig_vec, eig_val] = eig(L); 
    
    [eig_val, idx] = sort(diag(eig_val));
    eig_vec = eig_vec(:, idx);
    eig_vec = eig_vec(:, 2:end);
    
    eig_val = diag(eig_val(2:end));

end
