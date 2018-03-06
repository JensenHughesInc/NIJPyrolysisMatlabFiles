%% Running Script by Fengchang Yang
% 08/22/2017 Version 1.1
% Update from old version:
% --Enable pyrolysis module
% --Separate material properties file "MaterialProperties.m"
% --Boundary condition type and parameters are moved to running script
% --Mesh control parameters in script
% Polynomial Conductivity and Polynomial Heat Capacity

clear all;
close all;

lin_cond        = [2.9997e-4 4.405e-8 0.0];
poly_cond       = [2.9997e-4 4.4050e-8 1.4364e-25];
lin_heat_cap    = [1.0677 4.520e-5 0.0];
poly_heat_cap   = [1.0677 4.5200e-5 4.0694e-021];

%% Define Control Parameters
dt    = 0.5;      % time step (s)

filename    = 'poly_poly.mat';
xback       = 0.037;                        % total thickness of the material
nRegion     = 3;                            % total number of mesh refine region
samback1    = 0.006;                        % mesh refine region #1
samback2    = 0.012;                        % mesh refine region #2
mesh        = [0 samback1 samback2 xback];  
mesh_ne     = [20 10 10];                % element amount for each region

load(filename);

tend        = 200.;                         % end time of simulation
chcTisurf   = 300.;
chcTiback   = 300.;
Data = Tpro;

SA_flag = 'Nothing';

%% Define Boundary Conditions
% Handler for BCs
% Current supporting types:
% 'Flux' -- general flux involving both appled, convection, and radiation;
% 'Adiabatic' -- only heat flux on one side;
% 'Dirichlet' -- fixed wall temperature on one or each sides (implemented
% later;
BCtype = 'Flux';

% Define BCs according to handler
if(strcmp(BCtype, 'Flux'))
    q1vector = [0 25; tend 25];     % Applied heat flux at node 1 [kW]
    qnvector = [0 0; tend 0];       % Applied heat flux at node n [kW]
    hc1 = 0.01;                     % Heat transfer coefficient for convection at node 1 [kW/m^2/K]
    hcn = 0.01;                     % Heat transfer coefficient for convection at node n [kW/m^2/K]
    chce1 = 0.0;                    % Emissivity of both surfaces
elseif(strcmp(BCtype, 'Adiabatic'))
    q1vector = [0 0; tend 0];
    qnvector = [0 0; tend 0];
    hc1 = 0.0;
    hcn = 0.0;
    chce1 = 1.0;
elseif(strcmp(BCtype, 'Dirichlet'))
    
end

%% Run 1D FEM Simulation with Inputs
[curr_tdim, curr_tpro] = FEDM(dt, tend, xback, samback1, samback2, mesh_ne, chce1,...
    chcTisurf, chcTiback, q1vector, qnvector, hc1, hcn, SA_flag);

%% Post-processing
for i=1:length(Tpro(:,1))
    T1(i,1) = tdim(i);
    T1(i,2) = Tpro(i,1);
    
    T2(i,1) = tdim(i);
    T2(i,2) = Tpro(i,4);
end

for i=1:length(curr_tdim)
    tcurr = curr_tdim(i);
    currT1_data(i) = Decomp1D_InterpolateValue(T1, tcurr);
    currT2_data(i) = Decomp1D_InterpolateValue(T2, tcurr);
end

index = 1;
for i=1:length(curr_tdim)
    if mod(curr_tdim(i),10) == 0
        tdim_plot(index) = curr_tdim(i);
        ftemp_plot(index) = curr_tpro(i,1);
        mtemp_plot(index) = curr_tpro(i,3);
        btemp_plot(index) = curr_tpro(i,4);
        index = index + 1;
    end
end

figure;
hold on;
plot(curr_tdim, currT1_data-273.15, ':k');
plot(curr_tdim, currT2_data-273.15, '-.k');
plot(tdim_plot, ftemp_plot-273.15, '^k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
plot(tdim_plot, btemp_plot-273.15, 'vk', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
plot(tdim_plot, mtemp_plot-273.15, 'ok', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
H = legend('Frontside (Data)', 'Backside (Data)', 'Frontside (Predicted)',...
    'Backside (Predicted)','Middle (Predicted)', 'location', 'northwest');
xlabel('Time (sec)', 'FontSize', 12);
ylabel('Temperature (\circC)', 'FontSize', 12);
hold off;