function [] = Decomp1D_PlotTempProfilevTime( xdim_plot, Tall )
%   Plots the temperature profile versus time for all nodes
%
%    Calling Sequence:
%    [] = Decomp1D_PlotTempvTime( xdim_plot, Tall )
%
%    Input:   xdim_plot - locations of all nodes
%             Tall      - temperatures at all nodes at all time
%
%    Output:  figure
%

figure
plot(xdim_plot,Tall)
title('Temperatures Plotted over Time and Element Length')
xlabel('Distance from Surface (m)')
ylabel('Temperature (K)')