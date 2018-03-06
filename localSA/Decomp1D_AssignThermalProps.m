function [prop] = Decomp1D_AssignThermalProps( prop, ne, Teavg, node_div, aconst, bconst, cconst )
%   Assign thermal properties to a matrix based upon the material assigned
%   at each element. Another subroutine exists for the enthalpies.
%
%    Calling Sequence:
%    [prop] = Decomp1D_AssignThermalProps( prop, ne, Teavg, node_div, aconst, bconst, cconst )
%
%    Input:   prop   - initialized property from main code
%             ne     - number of elements in the system
%             Teavg  - number of materials in the system
%             node_div - the nodes located at the material divisions
%             aconst - the 'a' constant for calculation
%             bconst - the 'b' constant for calculation
%             cconst - the 'c' constant for calculation
%
%    Output:  prop - property matrix calculated using constants
%

for i = 1:ne
    for j = 1:length(node_div)
        if(i<node_div(j))
            % Initialize the thermal properties of the material as a function of
            % temperature (Teavg -> initial elemental temperatures)
            prop(i)  = (Teavg(i).^2*cconst(j)+Teavg(i)*bconst(j))+aconst(j);
            break
        end
    end
end