function returnValue = Decomp1D_InterpolateValue_2(vector, tcurr)


	returnValue = interp1qr(vector(:,1),vector(:,2),tcurr); % Custom Function for quicker interpolation
                                                            % 2013, Jose M. Mier
                                                            % Previously used interp1q, the
                                                            % custom function is a drop in replacement
    if isnan(returnValue)                                            
        if tcurr < vector(1,1)
            returnValue = vector(1,2);
        elseif tcurr > vector(length(vector(:,1)),1)
            returnValue = vector(length(vector(:,1)), 2);
        end
    end
                                                        
    
  % for i=1:length(vector(:,1))
  %     if tcurr == vector(i,1)
  %         returnValue = vector(i,2);
  %     elseif tcurr > vector(i,1) && tcurr < vector(i+1, 1)
  %         returnValue = (vector(i+1,2) - vector(i,2))/(vector(i+1,1) - vector(i,1))*(tcurr - vector(i,1)) + vector(i,2);
  %     end
  %  end
end
