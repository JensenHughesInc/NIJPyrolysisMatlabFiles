function [A] = Decomp1D_AssembleStiffnessMatrix( A, k, h1, hn, dx, ne )
%   Assemble the stiffness matrix for the system.
%
%    Calling Sequence:
%    [A] = Decomp1D_AssembleStiffnessMatrix( A, k, h1, hn, dx, ne )
%
%    Input:   A  - initialized stiffness matrix from main code
%             k  - the conduction coefficient matrix
%             h1 - convection coefficient for node 1
%             hn - convection coefficient for node n
%             dx - node spacing matrix
%             ne - number of elements in the system
%
%    Output:  A - assembled stiffness matrix
%
  for i = 1:ne
    A(i,i)      = A(i,i)+k(i)/dx(i);
    A(i,i+1)    = A(i,i+1)-k(i)/dx(i);
    A(i+1,i)    = A(i+1,i)-k(i)/dx(i);
    A(i+1,i+1)  = A(i+1,i+1)+k(i)/dx(i);
    % If this is the first node, add the heat transfer coefficient 
    % associated with this node
    if (i==1)
      A(i,i)    = A(i,i)+h1;
    end
    % If this is the last node, add the heat transfer coefficient
    % associated with this node
    if (i==ne)
      A(i+1,i+1)= A(i+1,i+1)+hn;
    end    
  end