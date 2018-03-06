function [mf] = Decomp1D_CalculateMassFlux_2( mf, xdim, drdt, n, ne, num_mtl )
%   Calculate the mass flux using the trapezoidal integration technique
%
%   S. Kraft Edited Version
%   June 14, 2016    
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
% for i = 1:(ne+(num_mtl-1))
%    xint    = xdim((n+(num_mtl-1)-i):n+(num_mtl-1)); % determines x-values from L to x'
%    drdtint = drdt((n+(num_mtl-1)-i):n+(num_mtl-1)); % determines dr/dt values from L to x'
%    mf(n+(num_mtl-1)-i) = trapz(xint,drdtint);       % computes the trapezoidal integral of 
%                                                         % drdtint WRT xint
% end
mf = flipud(-cumtrapz(xdim(end:-1:1),drdt(end:-1:1))); %This command has same output as loop above
