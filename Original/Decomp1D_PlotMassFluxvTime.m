function [] = Decomp1D_PlotMassFluxvTime( tdim, mfsurf )
%   Plots mass flux at the surface versus time
%
%    Calling Sequence:
%    [] = Decomp1D_MassFluxvTime( tdim, mfsurf )
%
%    Input:   tdim   - time matrix
%             mfsurf - mass flux at the surface for each time step
%
%    Output:  figure
%

figure
plot(tdim,mfsurf)
title('Mass Loss at Exposed Surface')
xlabel('Time (sec)')
ylabel('Mass Loss per unit Area at the Exposed Surface (kg/s-m^2)')