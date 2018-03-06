function [EDT] = Decomp1D_AssembleEDTMatrix( EDT, dh, dhg, dx, node_div, ne )
%   Assemble the transient mass matrix for the system.
%
%    Calling Sequence:
%    [EDT] = Decomp1D_AssembleEDTMatrix( EDT, dh, dhg, dx, node_div, ne )
%
%    Input:   EDT  - initialized density v time matrix from main code
%             dh   - the enthalpy matrix
%             dhg  - the enthalpy of gas matrix
%             dx - node spacing
%             node_div - the nodes located at material intersections
%             ne - number of elements in the system
%
%    Output:  EDT - assembled density versus time matrix
%

  nctr = 0;
  for i = 1:ne
      if(i==node_div(nctr+1))
          nctr = nctr+1;
      end
      EDT(i+nctr,i+nctr)     = EDT(i+nctr,i+nctr)+(dh(i)-dhg(i))*dx(i)/3;
      EDT(i+nctr,i+1+nctr)   = EDT(i+nctr,i+1+nctr)+(dh(i)-dhg(i))*dx(i)/6;
      EDT(i+1+nctr,i+nctr)   = EDT(i+1+nctr,i+nctr)+(dh(i)-dhg(i))*dx(i)/6;
      EDT(i+1+nctr,i+1+nctr) = EDT(i+1+nctr,i+1+nctr)+(dh(i)-dhg(i))*dx(i)/3;
  end