function [Feavg] = Decomp1D_CalculateFeavg( Fea1, Fea2, node_div, n )
%   Generate the x-dimension matrix for future calculations
%
%    Calling Sequence:
%    [xdim] = Decomp1D_GenerateXDim( xdim, num_mtl, node_div, n, dx )
%
%    Input:   Fea1     - Feavg for node 1 of each element
%             Fea2     - Feavg for node 2 of each element
%             node_div - the nodes located at the material divisions
%             n        - number of nodes in the system
%
%    Output:  Feavg - average mass fraction for each element
%
  nctr = 0.0;
  for i = 1.0:n
      if(i==node_div(nctr+1.0))
          if(node_div(nctr+1.0)==n)
              break
          end
          Fea1(i-nctr) = [];     % eliminate F for material i
          Fea2(i-nctr) = [];     % eliminate F for material j
          nctr = nctr+1.0;
      end
  end
  Feavg = (Fea1+Fea2)/2;        % calculate the elemental avg fractions