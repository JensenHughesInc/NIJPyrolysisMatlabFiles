function [Ffrac] = Decomp1D_UpdateMassFraction( Ffrac, r, rhod, drho, n, node_div,count )
%   Calculate the instantaneous mass fractions at the nodes based upon the
%   respective change in density at each node
%
%    Calling Sequence:
%    [Ffrac] = Decomp1D_UpdateMassFraction( Ffrac, r, rhod, drho, node_div )
%
%    Input:   Ffrac  - initialized mass fraction matrix from main code
%             r      - densities from the current time step
%             rhod   - the decomposed density
%             drho   - (rhov-rhod) - the difference of virgin and decomposed densities
%             n      - number of nodes
%             node_div - the intersection location of the materials
%
%    Output:  Ffrac - instantaneous mass fraction matrix
%
for i = 1:n+(length(node_div)-1)
    for j = 1:length(node_div)
        if(i<=node_div(j) && drho(j) > 0)
            % Calculate the instantaneous mass fraction
            Ffrac(i) = (r(i)-rhod(j))/drho(j);
            break
        end
    end
end