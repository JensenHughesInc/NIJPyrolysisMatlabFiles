function [xdim] = Decomp1D_GenerateXDim( xdim, num_mtl, node_div, n, dx )
%   Generate the x-dimension matrix for future calculations
%
%    Calling Sequence:
%    [xdim] = Decomp1D_GenerateXDim( xdim, num_mtl, node_div, n, dx )
%
%    Input:   xdim     - initialized x-dims from main code
%             num_mtl  - number of materials in the system
%             node_div - the nodes located at the material divisions
%             n        - number of nodes in the system
%             dx       - matrix of individual element lengths
%
%    Output:  xdim - x-dimension matrix
%

nctr = 0;
dxctr = 1;
for i=1:n+(num_mtl-1)
    if(i==1)
        xdim(i) = 0;
    elseif((nctr~=0)&&(i-1==node_div(nctr)))
        xdim(i) = xdim(i-2)+dx(dxctr);
        dxctr = dxctr+1;
    else   
        xdim(i) = xdim(i-1)+dx(dxctr);
        if(i==node_div(nctr+1))
            nctr = nctr+1;
            if(nctr==length(node_div))
                nctr = nctr-1;
            end
        else
            dxctr = dxctr+1;
        end
    end
end
    