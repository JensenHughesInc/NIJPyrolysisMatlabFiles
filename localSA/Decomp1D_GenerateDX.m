function [dx] = Decomp1D_GenerateDX( dx, mesh_ne, mesh_dx, ne )
%   Generate the node spacing matrix
%
%    Calling Sequence:
%    [dx] = Decomp1D_GenerateDX( dx, mesh_ne, mesh_dx )
%
%    Input:   dx      - initialized dx from main code
%             mesh_ne - number of elements per refinement zone
%             mesh_dx - the node spacing in each refinement zone
%
%    Output:  dx - node spacing matrix
%

dxctr   = 1;                % initialize counter for tracking location in [mesh_dx]
curr_ne = mesh_ne(dxctr);   % initialize current number of elements to track to
for i = 1:ne
    dx(i) = mesh_dx(dxctr);
    if(i==curr_ne)
        dxctr = dxctr + 1;
        if(dxctr<=length(mesh_ne))
            curr_ne = curr_ne + mesh_ne(dxctr);
        end
    end
end