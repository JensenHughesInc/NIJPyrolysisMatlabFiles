function [Tinit] = Decomp1D_InterpolateTinit( Tinitvector, dx, n )
%   Calculate the initial temperatures throughout the sample using linear
%   interpolation
%
%    Calling Sequence:
%    [Tinit] = Decomp1D_InterpolateTinit( Tinitvector, xdim )
%
%    Input:   Tinitvector - interpolation vector for Tinit
%             xdim        - x-dimension matrix
%
%    Output:  Tinit - initial temperature matrix
%
Tinit = zeros(n,1);

x_loc = 0.0;
for i = 1:n
    if i == 1
        Tinit(i) = Tinitvector(1,2);
    elseif i == n
        Tinit(i) = Tinitvector(length(Tinitvector),2);
    else
        for j = 1:length(Tinitvector)
            if(x_loc<=Tinitvector(j,1))
                Tinit(i) = Tinitvector(j,2)-((Tinitvector(j,2)-Tinitvector(j-1,2))/...
                    (Tinitvector(j,1)-Tinitvector(j-1,1)))*(Tinitvector(j,1)-x_loc);
    
                break
            end
        end
    end
    if(i~=n)
        x_loc = x_loc + dx(i);
    end
end


