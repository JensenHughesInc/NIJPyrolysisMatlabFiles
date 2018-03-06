function [EDT_drdt] = Decomp1D_PremultiplyEDTdrdt( EDT_drdt, EDT, drdt, node_div, n )
%   Premultiply the EDT and drdt matrices for multiple materials
%
%    Calling Sequence:
%    [EDT_drdt] = Decomp1D_PremultiplyEDTdrdt( EDT_drdt, EDT, drdt, node_div, ne )
%
%    Input:   EDT_drdt - matrix initialized in main code
%             EDT  - initialized density v time matrix from main code
%             drdt - the change in density with time matrix
%             node_div - the nodes located at material intersections
%             ne - number of elements in the system
%
%    Output:  EDT_drdt - premultiplied EDT and drdt matrix
%

% Calculate the current EDT and drdt matrix multiplication
premult = EDT*drdt;

% Find material divisions in premutliplied matrix and sum together,
% reducing overall matrix length from (n+num_mtl-1) to n
nctr = 0;
for i = 1:n
    EDT_drdt(i) = EDT_drdt(i)+premult(i+nctr);
    if(i==node_div(nctr+1))
        nctr = nctr+1;
        if(i~=node_div(length(node_div)))
            EDT_drdt(i) = EDT_drdt(i)+premult(i+nctr);
        end
    end
end