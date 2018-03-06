function [] = Decomp1D_PlotMassFluxvTime_ValidationData( filename, tdim_range, mpro_range, tdim, mfsurf, ymin, ymax )
%   Plots mass flux at the surface versus time
%
%    Calling Sequence:
%    [] = Decomp1D_PlotMassFluxvTime_ValidationData( filename, tdim_range, mpro_range, tdim, mfsurf, ymin, ymax )
%
%    Input:   filename   - name of .xls file containing validation data
%             tdim_range - range of cells in file containing time data
%             mpro_range - range of cells in file containing surface mass flux data
%             tdim       - time step data from model
%             mfsurf     - mass flux at the surface for each time step
%             ymin       - minimum y-axis value for the figure
%             ymax       - maximum y-axis value for the figure
%
%    Output:  figure
%

% Extract data from stored validation data files
val_tdim    = xlsread(filename,1,tdim_range);
val_mfsurf  = xlsread(filename,1,mpro_range);

figure
% Plot the surface mass flux v time from code
plot(tdim,mfsurf)
hold
% % Plot the surface mass flux v time from Validation data
plot(val_tdim,val_mfsurf(:,1),'r');
axis([0 2000 ymin ymax])
legend('FE Model','Validation Data',4)
title('Mass Flux at Surface v Time')
xlabel('Time (sec)')
ylabel('Mass Flux per unit Area at the Exposed Surface (kg/s-m^2)')