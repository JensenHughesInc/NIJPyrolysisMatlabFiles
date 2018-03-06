function [B] = Decomp1D_AssembleTransientMatrix( B, r, c, dx, node_div, ne )
%   Assemble the transient mass matrix for the system.
%
%    Calling Sequence:
%    [B] = Decomp1D_AssembleTransientMatrix( B, r, c, dx, node_div, ne )
%
%    Input:   B  - initialized transient mass matrix from main code
%             r  - the density matrix
%             c  - specific heat matrix
%             dx - node spacing
%             node_div - the nodes located at material intersections
%             ne - number of elements in the system
%
%    Output:  B - assembled transient mass matrix
%


  nctr = 0;
  for i = 1:ne
      if(i==node_div(nctr+1))
          nctr = nctr+1;
      end
      B(i,i)     = B(i,i)+(3*r(i+nctr)+r(i+1+nctr))*c(i)*dx(i)/12;
      B(i,i+1)   = B(i,i+1)+(r(i+nctr)+r(i+1+nctr))*c(i)*dx(i)/12;
      B(i+1,i)   = B(i+1,i)+(r(i+nctr)+r(i+1+nctr))*c(i)*dx(i)/12;
      B(i+1,i+1) = B(i+1,i+1)+(r(i+nctr)+3*r(i+1+nctr))*c(i)*dx(i)/12;
  end
  