function H = heat_kernel(A,V,E,t)
   %[V, E] = laplacian(A);
   H = V*diag(diag(exp(-E*scale(t))))*V';
end