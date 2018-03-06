function [mf] = Decomp1D_CalculateMassFlux( mf, xdim, drdt, n, ne, num_mtl )
%   Calculate the mass flux using the trapezoidal integration technique
%
%    Calling Sequence:
%    [mf] = Decomp1D_CalculateMassFlux( mf, xdim, drdt, n, ne, num_mtl )
%
%    Input:   mf   - initialized mass flux from main code
%             xdim - x-dimension matrix associated with densities
%             drdt - change of density with time of current time step
%             n    - number of nodes
%             ne   - number of elements
%             num_mtl  - number of materials
%
%    Output:  mf - mass flux at the nodes
%

% Mass flux calculation
for i = 1:(ne+(num_mtl-1))
   xint    = xdim((n+(num_mtl-1)-i):n+(num_mtl-1)); % determines x-values from L to x'
   drdtint = drdt((n+(num_mtl-1)-i):n+(num_mtl-1)); % determines dr/dt values from L to x'
%    if(length(xint)>100)
%        mf(n+(num_mtl-1)-i) = gqint(xint,drdtint);
%    else
%        mf(n+(num_mtl-1)-i) = trapz(xint,drdtint);
%    end       
    mf(n+(num_mtl-1)-i) = trapz(xint,drdtint);
end

