function [CONV] = Decomp1D_AssembleConvectionMatrix( CONV, mf, cg, node_div, ne )
%   Assemble the transient mass matrix for the system.
%
%    Calling Sequence:
%    [CONV] = Decomp1D_AssembleConvectionMatrix( CONV, mf, cg, node_div, ne)
%
%    Input:   CONV     - initialized convection matrix from main code
%             mf       - the mass flux matrix
%             cg       - specific heat of gas matrix
%             node_div - the nodes located at material intersections
%             ne       - number of elements in the system
%
%    Output:  CONV - assembled convenction matrix
%

  nctr = 0;
  for i = 1:ne
      if(i==node_div(nctr+1))
          nctr = nctr+1;
      end
      CONV(i,i)       = CONV(i,i)+(-mf(i+nctr)/3-mf(i+1+nctr)/6)*cg(i);
      CONV(i,i+1)     = CONV(i,i+1)+(mf(i+nctr)/3+mf(i+1+nctr)/6)*cg(i);
      CONV(i+1,i)     = CONV(i+1,i)+(-mf(i+nctr)/6-mf(i+1+nctr)/3)*cg(i);
      CONV(i+1,i+1)   = CONV(i+1,i+1)+(mf(i+nctr)/6+mf(i+1+nctr)/3)*cg(i);
  end 