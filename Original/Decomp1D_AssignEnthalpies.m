function [prop] = Decomp1D_AssignEnthalpies( prop, ne, Teavg, Tr, node_div, aconst, bconst, cconst, dhrec, trans, Feavg )
%   Assign enthalpies to a matrix based upon the material assigned
%   at each element. Another subroutine exists for the thermal properties.
%
%    Calling Sequence:
%    [prop] = Decomp1D_AssignEnthalpies( prop, ne, Teavg, node_div, aconst, bconst, cconst )
%
%    Input:   prop   - initialized property from main code
%             ne     - number of elements in the system
%             Teavg  - number of materials in the system
%             Tr     - reference temperature
%             node_div - the nodes located at the material divisions
%             aconst - the 'a' constant for calculation
%             bconst - the 'b' constant for calculation
%             cconst - the 'c' constant for calculation
%             dhrec  - heat of reaction
%             trans  - transition points (F) for decompostion phases
%             Feavg  - avg instantaneous mass fraction of element
%
%    Output:  prop - property matrix calculated using constants
%

pctr = 1;   % counter to track the current phase
for i = 1:ne
    for j = 1:length(node_div)
        if(i<node_div(j))
            % Check to determine current decomposition stage
            for k = 1:size(trans,2)
                if(Feavg(i)>0 && Feavg(i)<trans(j,k))
                    pctr = pctr + 1;
                end
            end
            % Initialize the enthalpy (h) and enthalpy of gas (hg)
            %   - Properties for dh are based on virgin material b/c enthalpy is
            %       assumed to not be a function of decomposition
            prop(i)  = ((Teavg(i).^3-Tr(i).^3)*cconst(j)/3+(Teavg(i).^2-...
                Tr(i).^2)*bconst(j)/2+(Teavg(i)-Tr(i))*aconst(j))+dhrec(j,pctr);
            pctr = 1;   % reset decomposition phase counter
            break
        end
    end
end