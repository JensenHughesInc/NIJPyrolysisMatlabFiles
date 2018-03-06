function [r] = Decomp1D_InitializeDensity( r, n, num_mtl, node_div, rhov )
%   Assign densities in the density matrix according to appropriate
%   materials and locations
%
%    Calling Sequence:
%    [r] = Decomp1D_InitializeDensity( num_mtl, node_div, rhov )
%
%    Input:   r        - initialized density matrix from main code
%             n        - the number of nodes in the system
%             num_mtl  - number of materials in the system
%             node_div - the nodes at the divisions of the materials
%             rhov     - the virgin densities of all materials
%
%             NOTE: The first node division will always be zero
%
%    Output:  r - density matrix with dims (n+1,1)
%

for i = 1:n+(num_mtl-1)
    for j = 1:length(node_div)
        if(i<=node_div(j))
            r(i) = rhov(j);     % virgin material density matrix
            break
        end
        if(i>=node_div(length(node_div)))
            r(i) = rhov(j);     % define final density in last material
        end
    end
end