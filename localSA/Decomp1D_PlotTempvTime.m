function [] = Decomp1D_PlotTempvTime( tdim, Tpro, legendFlag )
%   Plots the temperatures versus time for the model
%
%    Calling Sequence:
%    [] = Decomp1D_PlotTempvTime( tdim, Tpro )
%
%    Input:   tdim - time matrix
%             Tpro - temperatures recorded over time at selected nodes
%
%    Output:  figure
%

figure
plot(tdim,Tpro)
if legendFlag == 'A'
    legend('Exposed Surface', 'Mid Laminate', 'Unexposed Surface', 'Location', 'NorthEastOutside')
elseif legendFlag == 'B'
    legend('Insulation Surface', 'Laminate Surface', 'Mid Laminate', 'Unexposed Surface', 'Location', 'NorthEastOutside');
elseif legendFlag == 'C'
    legend('Exposed Surface', 'Interface One', 'Center Balsa Wood', 'Interface Two', 'Unexposed Surface');
elseif legendFlag == 'D'
    legend('Exposed Surface', 'Mid Laminate', 'Laminate/Superwool', 'Mid Superwool', 'Unexposed Surface');
elseif legendFlag == 'E'
    legend('Exposed Surface', 'Superwool/Laminate', 'Mid Laminate', 'Laminate/Superwool', 'Unexposed Surface');
elseif legendFlag == 'F'
    legend('Exposed Surface', 'Laminate/Balsa Wood', 'Mid Balsa Wood', 'Balsa Wood/Laminate', 'Laminate/Superwool', 'Unexposed Surface');
end
title('Temperature Profiles v Time')
xlabel('Time (sec)')
ylabel('Temperature (K)')
