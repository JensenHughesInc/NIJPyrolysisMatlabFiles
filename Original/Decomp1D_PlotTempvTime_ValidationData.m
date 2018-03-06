function [] = Decomp1D_PlotTempvTime_ValidationData( filename, tdim_range, Tpro_range, tdim, Tpro )
%   Plots the temperatures versus time for the model and Validation data
%
%    Calling Sequence:
%    [] = Decomp1D_PlotMassFluxvTime_ValidationData( filename, tdim_range, Tpro_range, tdim, Tpro )
%
%    Input:   filename   - name of .xls file containing validation data
%             tdim_range - range of cells in file containing time data
%             Tpro_range - range of cells in file containing temperature data
%             tdim       - time step data from model
%             Tpro       - temperatures at selected nodes over time
%
%    Output:  figure
%

% Extract data from stored validation data files
val_tdim    = xlsread(filename,1,tdim_range);
val_Tpro    = xlsread(filename,1,Tpro_range);

figure
% Plot the T profile v time from code
plot(tdim,Tpro)
hold
% % Plot the T profile v time from Validation Data
plot(val_tdim,val_Tpro(:,1),'-.b',val_tdim,val_Tpro(:,2),'-.g',...
    val_tdim,val_Tpro(:,3),'-.r');
legend('FE Model (x=0mm)','FE Model (x=6.35mm)', 'FE Model (x=12.7mm)',...
    'Validation (x=0mm)','Validation (x=6.35mm)','Validation (x=12.7mm)',4)
title('Temperature Profiles v Time')
xlabel('Time (sec)')
ylabel('Temperature (K)')