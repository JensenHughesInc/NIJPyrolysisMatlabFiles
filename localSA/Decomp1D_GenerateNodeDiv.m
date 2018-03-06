function [node_div] = Decomp1D_GenerateNodeDiv( node_div, num_mtl, mtl_div, mesh, mesh_ne )
%   Generate the material division node matrix
%
%    Calling Sequence:
%    [node_div] = Decomp1D_GenerateNodeDiv( node_div, num_mtl, mtl_div, mesh, mesh_ne )
%
%    Input:   node_div - initialized node_div matrix from main code
%             num_mtl  - number of materials in the model
%             mtl_div  - the Cartesian coordinate locations of divisions
%             mesh     - the local mesh 'refinement' divisions
%             mesh_ne  - the number of elements in each mesh refinement 'zone'
%
%    Output:  node_div - nodes located at material divisions
%

numelem = 0;
for i = 1:num_mtl
    curr_div = mtl_div(i);
    for j = 1:length(mesh)
        if(curr_div==mesh(j))
            for k = 1:(j-1)
                numelem = numelem + mesh_ne(k);
            end
            node_div(i) = numelem + 1;
            numelem = 0;
        end
    end
end