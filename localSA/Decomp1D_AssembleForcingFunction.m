function [F] = Decomp1D_AssembleForcingFunction( F, Q, q1, qn, h1, hn, To1, Ton, dx, ne )
%   Assemble the transient mass matrix for the system.
%
%    Calling Sequence:
%    [F] = Decomp1D_AssembleForcingFunction( F, Q, q1, qn, h1, hn, To1,
%    Ton, dx, ne )
%
%    Input:   F   - initialized forcing function from main code
%             Q   - internal heat generation
%             q1  - heat flux at node 1
%             qn  - heat flux at node n
%             h1  - convection coefficient at node 1
%             hn  - convection coefficient at node n
%             To1 - gas temperature at node 1
%             Ton - gas temperature at node n
%             dx  - node spacing
%             ne  - number of elements in the system
%
%    Output:  F - assembled forcing function matrix
%

  for i = 1:ne 
    F(i)    = F(i)+dx(i)*Q/2;      % heat generation term
    F(i+1)  = F(i+1)+dx(i)*Q/2;    % heat generation term
    if (i==1)
      F(i)  = F(i)+q1+h1*To1;   % heat flux at node 1
    end
    if (i==ne)
      F(i+1)= F(i+1)+qn+hn*Ton; % heat flux at node n
    end
  end 