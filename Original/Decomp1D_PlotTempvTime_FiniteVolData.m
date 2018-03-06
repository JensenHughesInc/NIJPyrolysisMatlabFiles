function [] = Decomp1D_PlotTempvTime_FiniteVolData( filename, tdim_range, Tpro_range, tdim, Tpro )
%   Plots temperature profile versus time for model and Finite Vol data
%
%    Calling Sequence:
%    [] = Decomp1D_PlotMassFluxvTime_FiniteVolData( filename, tdim_range, Tpro_range, tdim, Tpro )
%
%    Input:   filename   - name of .xls file containing validation data
%             tdim_range - range of cells in file containing time data
%             Tpro_range - range of cells in file containing temperature data
%             tdim       - time step data from model
%             mfsurf     - mass flux at the surface for each time step
%
%    Output:  figure
%

% Extract data from stored validation data files
fv_tdim     = xlsread(filename,1,tdim_range);
fv_Tpro     = xlsread(filename,1,Tpro_range);

figure
% Plot the T profile v time from code
plot(tdim,Tpro)
hold
% Plot the T profile v time from the Finite Vol model
plot(fv_tdim,fv_Tpro(:,1),'-.b',fv_tdim,fv_Tpro(:,2),'-.g',...
    fv_tdim,fv_Tpro(:,3),'-.r');
legend('FE - Node(1)','FE - Node(n/2)','FE - Node(n)',...
    'FV - Node(1)','FV - Node(n/2)','FV - Node(n)',4)
title('Temperature Profiles v Time for Composite')
xlabel('Time (sec)')
ylabel('Temperature (K)')