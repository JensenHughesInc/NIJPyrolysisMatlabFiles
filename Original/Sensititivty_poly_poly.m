%% Sensitivity Analysis File for poly_poly case
% --Fengchang Yang, 08/22/2017
%-------------------------------------------------------------------------%
%Discription

clear all;
close all;
%% Parameters define
Para = {'rhov1', 'rhod1', 'E', 'n'};

%% Define Control Parameters
dt    = 0.5;      % time step (s)

filename    = 'poly_poly.mat';

load(filename);

xback       = 0.037;                        % total thickness of the material
nRegion     = 3;                            % total number of mesh refine region
samback1    = 0.006;                        % mesh refine region #1
samback2    = 0.012;                        % mesh refine region #2
mesh        = [0 samback1 samback2 xback];  
mesh_ne     = [20 10 10];                   % element amount for each region

tend        = 1000.;                        % end time of simulation
chcTisurf   = 300.;
chcTiback   = 300.;

SA = 1.01;

T = {};

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
    q1vector = [0 25; tend 25];     % Applied heat flux at node 1 [kW/m^2]
    qnvector = [0 0; tend 0];       % Applied heat flux at node n [kW/m^2]
    hc1 = 0.01;                     % Heat transfer coefficient for convection at node 1 [kW/m^2/K]
    hcn = 0.01;                     % Heat transfer coefficient for convection at node n [kW/m^2/K]
    chce1 = 0.0;                    % Emissivity of both surfaces
elseif(strcmp(BCtype, 'Adiabatic'))
    q1vector = [0 0; tend 0];
    qnvector = [0 0; tend 0];
    hc1 = 0.0;
    hcn = 0.0;
    chce1 = 0.0;
elseif(strcmp(BCtype, 'Dirichlet'))
    
end
%% Standard run

SA_flag = 'Standard';
[curr_tdim, curr_T] = FEDM(dt, tend, xback, samback1, samback2, mesh_ne, chce1,...
    chcTisurf, chcTiback, q1vector, qnvector, hc1, hcn, SA_flag, SA);
T{1} = curr_T;

%% Local SA run
for i = 1:length(Para)
    SA_flag = Para{i};
    [curr_tdim, curr_T] = FEDM(dt, tend, xback, samback1, samback2, mesh_ne, chce1,...
        chcTisurf, chcTiback, q1vector, qnvector, hc1, hcn, SA_flag, SA);
    SA_flag
    T{i+1} = curr_T;
end

%% Post-processing
figure;
hold on;
plot(curr_tdim, (T{2}(:,1)-T{1}(:,1))/(SA-1)./(T{1}(:,1)), '-', 'linewidth', 1.5);
plot(curr_tdim, (T{3}(:,1)-T{1}(:,1))/(SA-1)./(T{1}(:,1)), '--', 'linewidth', 1.5);
plot(curr_tdim, (T{4}(:,1)-T{1}(:,1))/(SA-1)./(T{1}(:,1)), '-.', 'linewidth', 1.5);
plot(curr_tdim, (T{5}(:,1)-T{1}(:,1))/(SA-1)./(T{1}(:,1)), ':', 'linewidth', 1.5);
H = legend(Para{1}, Para{2}, Para{3}, Para{4}, 'location', 'southwest');
xlabel('Time (sec)', 'FontSize', 12);
ylabel('Normalized SA coeff', 'FontSize', 12);
hold off;

figure;
hold on;
plot(curr_tdim, T{2}(:,1), '-', 'linewidth', 1.5);
plot(curr_tdim, T{3}(:,1), '--', 'linewidth', 1.5);
plot(curr_tdim, T{4}(:,1), '-.', 'linewidth', 1.5);
plot(curr_tdim, T{5}(:,1), ':', 'linewidth', 1.5);
plot(tdim(1:10:end), Tpro(1:10:end,1), '^', 'MarkerEdgeColor', 'k');
H = legend(Para{1}, Para{2}, Para{3}, Para{4}, 'Exp', 'location', 'southeast');
xlabel('Time (sec)', 'FontSize', 12);
ylabel('Temperature (K)', 'FontSize', 12);
hold off;
