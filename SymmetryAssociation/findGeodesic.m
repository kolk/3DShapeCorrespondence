function geoDist = findGeodesic(vertices, shape)
    size(vertices)
    geoDist = zeros(length(vertices), 1); 
    for cntr = 1:length(vertices)
        faces = shape.surface.TRIV;
        vertex_position = [shape.surface.X shape.surface.Y shape.surface.Z];
        start_vertex = vertices(cntr,1);
               
        %{
        %  D distances from id, # vertices by 1 matrix
        %  G exact gradients of D(in 3d), # vertices by 3 matrix
        [D,G] = exactgeodesic(vertex_position, faces, start_vertex);
        %}
        
       % perform the front propagation, Q contains an approximate segementation
       [D,S,Q] = perform_fast_marching_mesh(vertex_position, faces, start_vertex);
       geoDist(cntr) = D(vertices(cntr, 2));
    end
end