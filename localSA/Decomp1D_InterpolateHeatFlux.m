function [q1, qn] = Decomp1D_InterpolateHeatFlux( q1vector, qnvector, tcurr )
%   Calculate the heat flux at the current timestep using two interpolation
%   vectors
%
%    Calling Sequence:
%    [q1, qn] = Decomp1D_InterpolateHeatFlux( q1vector, qnvector, tcurr )
%
%    Input:   q1vector - interpolation vector for q1
%             qnvector - interpolation vector for qn
%             tcurr    - current time
%
%    Output:  q1 - heat flux at node 1
%             qn - heat flux at node n
%

for i = 1:length(q1vector)
    % Find location in interpolation vector of current time step
    if(tcurr<=q1vector(i,1))
        if(tcurr==0)
            q1 = q1vector(1,2);
        else
            q1 = q1vector(i,2)-((q1vector(i,2)-q1vector(i-1,2))/...
                (q1vector(i,1)-q1vector(i-1,1)))*(q1vector(i,1)-tcurr);
        end
        break
    end
end

for i = 1:length(qnvector)
    % Find location in interpolation vector of current time step
    if(tcurr<=qnvector(i,1))
        if(tcurr==0)
            qn = qnvector(1,2);
        else
            qn = qnvector(i,2)-((qnvector(i,2)-qnvector(i-1,2))/...
                (qnvector(i,1)-qnvector(i-1,1)))*(qnvector(i,1)-tcurr);
        end
        break
    end
end

