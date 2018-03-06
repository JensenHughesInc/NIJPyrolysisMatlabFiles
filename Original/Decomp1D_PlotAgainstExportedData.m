function [] = Decomp1D_PlotAgainstExportedData( filename, tdim_range, mpro_range, Tpro_range, tdim, mfsurf, Tpro )
%   Plots the current model data against previously exported model data
%
%    Calling Sequence:
%    [] = Decomp1D_PlotMassFluxvTime_ValidationData( filename, tdim_range, Tpro_range, tdim, Tpro )
%
%    Input:   filename   - name of .xls file containing validation data
%             tdim_range - range of cells in file containing time data
%             mpro_range - range of cells in file containing mass flux data
%             Tpro_range - range of cells in file containing temperature data
%             tdim       - time step data from model
%             mfsurf     - mass flux at the surface over time
%             Tpro       - temperatures at selected nodes over time
%
%    Output:  figure
%

% Extract data from stored validation data files
data_tdim   = xlsread(filename,1,tdim_range);
data_mfsurf = xlsread(filename,1,mpro_range);
data_Tpro   = xlsread(filename,1,Tpro_range);

figure
% Plot the current surface mass flux v time
plot(tdim,mfsurf)
hold
% Plot the surface mass flux v time from exported data set
plot(data_tdim,data_mfsurf(:,1),'r');
legend('Uniform Mesh','Non-uniform Mesh',4)
title('Mass Flux at Surface v Time')
xlabel('Time (sec)')
ylabel('Mass Flux per unit Area at the Exposed Surface (kg/s-m^2)')

figure
% Plot the T profile v time from code
plot(tdim,Tpro)
hold
% Plot the T profile v time from the Finite Vol model
plot(data_tdim,data_Tpro(:,1),'-.b',data_tdim,data_Tpro(:,2),'-.g',...
    data_tdim,data_Tpro(:,3),'-.r');
legend('Uniform Mesh - Node(1)','Uniform Mesh - Node(n/2)','Uniform Mesh - Node(n)',...
    'Non-uniform Mesh - Node(1)','Non-uniform Mesh - Node(n/2)','Non-uniform Mesh - Node(n)',4)
title('Temperature Profiles v Time for Composite')
xlabel('Time (sec)')
ylabel('Temperature (K)')