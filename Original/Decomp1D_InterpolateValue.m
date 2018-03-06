function returnValue = Decomp1D_InterpolateValue(vector, tcurr)

if tcurr < vector(1,1)
    returnValue = vector(1,2);
elseif tcurr > vector(length(vector(:,1)),1)
    returnValue = vector(length(vector(:,1)), 2);
else
	returnValue = interp1(vector(:,1),vector(:,2),tcurr);
    
  % for i=1:length(vector(:,1))
  %     if tcurr == vector(i,1)
  %         returnValue = vector(i,2);
  %     elseif tcurr > vector(i,1) && tcurr < vector(i+1, 1)
  %         returnValue = (vector(i+1,2) - vector(i,2))/(vector(i+1,1) - vector(i,1))*(tcurr - vector(i,1)) + vector(i,2);
  %     end
  %  end
end
