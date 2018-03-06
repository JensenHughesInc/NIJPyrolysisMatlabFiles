function [drdt] = Decomp1D_CalculateDRDT( drdt, r, drho, rhod, nak, Aak, Eak, trans, Ffrac, Rbar, T, inert, n, node_div)
%   Calculate the change in density based upon the last time step
%
%    Calling Sequence:
%    [drdt] = Decomp1D_CalculateDRDT( drdt, r, drho, rhod, nak, Aak, Eak, trans, Ffrac, Rbar, T, n, node_div, dt )
%
%    Input:   drdt - initiliazed density change matrix
%             r    - density matrix from last time step
%             drho - (rhov-rhod) - difference between virgin and decomposed
%             rhod - decomposed density matrix
%             nak  - density fraction power arrhenius kinetics
%             Aak  - pre-exponential factor arrhenius kinetics
%             Eak  - activation energy arrhenius kinetics
%             trans- transition points (F) for decomposition
%             Ffrac- instantaneous mass fraction at nodes
%             Rbar - ideal gas constant
%             T    - current nodal temperatures
%             inert- inert material matrix
%             n    - number of nodes
%             node_div - the nodes located at material intersections
%
%    Output:  drdt - change in density with time
%
 drdt_new = 0.0;    % variable to store newly calculated drdt in for comparison

 nctr = 1;  % counter to track material
 pctr = 1;  % counter to track current decomposition phase
 for i = 1:n+(length(node_div)-1)
     % Determine if pctr should be advanced to the next transition point
     for j = 1:size(trans,2)
         if(Ffrac(i)>0 && Ffrac(i)<trans(nctr,j))
             pctr = pctr+1;
         end         
     end
     
     if(inert(nctr)==0)
         % Calculate the change in density
         drdt_new = -1*drho(nctr)*(((r(i)-rhod(nctr))/drho(nctr))^nak(nctr,pctr))...
            *Aak(nctr,pctr)*exp(-Eak(nctr,pctr)/Rbar/T(i-nctr+1));
        % Determine if drdt is will result in a positive change in density
        if(drdt_new>0)
            drdt(i) = 0;
        else
            drdt(i) = drdt_new;
        end
        % Determine and correct if calculated change in density will cause
        % the nodal density to drop below the decomposed density
        if((r(i)+drdt(i))<rhod(nctr))
            drdt(i) = -1*(r(i)-rhod(nctr));
        end
     else
         drdt(i) = 0;
     end
     
     % Determine if nctr should be advanced to the next material
     if(i==node_div(nctr) && i~=node_div(length(node_div))) 
         nctr = nctr+1;         
     end
     pctr = 1;  % Reset the phase transition point counter
 end

