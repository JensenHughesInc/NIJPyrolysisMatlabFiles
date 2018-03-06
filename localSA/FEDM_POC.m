%{
    Functional version of the Finite Element Decomposition Model:
    Inputs:
        * input_file    -  The input file with the geometry and boundary condition data
        * dt            -  The time step used
        * tend          -  The final time for the scenario
%}
% function [tdim, Tpro] = FEDM_POC(akv1, bkv1, ckv1, acv1, bcv1, ccv1, dt, tend,...
%    xback, samback, chce1, chcrho, chch, chcTisurf, chcTiback, chcF, FileName)

function [tdim, Tpro] = FEDM_POC(akv1, bkv1, ckv1, acv1, bcv1, ccv1)

format short e

%%  Documentation and Versions

% VERSION 1.0
%   - Solves the general equation described above
%   - Linear elements developed using the Galerkin Formulation

% VERSION 1.1
%   - Boundary Condition issue corrected
%   - Matrix solution issue (matrices from previous time step) corrected
%   - Internal clock added to determine total elapsed time of code

% VERSION 1.2
%   - Decomposition terms added to the model, including the instantaneous
%       mass fraction, to track the material decomposition and respective
%       property changes
%   - For a constant volume element: F = (p - pd)/(pv - pd)
%   - Properties are dependent on decomposition as follow:
%           k = kv*F+(1-F)*kc
%           c = cv*F+(1-F)*cc
%   - NOTE: k & C is calculated using avg F over element

% VERSION 1.3
%   - Ability to model mutliple materials added to model
%       - densities are dually defined at material intersection
%   - Subroutines added to simplify code
%   - Allow use of varying external gas temperatures with time
%   - Heat flux can be varied with time using linear 

% VERSION 1.4
%   - Ability to locally refine mesh
%   - Can define multiple heats of decomposition to support multiple
%       materials
%   - Can define tranistion points for Arrhenius kinetics decomposition

%%  Initialize CPU Time Counter

time_init = cputime;    % initializes a time to determine the total
                            % elapsed time

%%  Initialize Time Values
% The following variables define the total time to run the analysis and the
% time step to use during analysis

dt    = 0.5;       % time step (s)
tend  = 200.;     % end time for simulation (s)
nt    = tend/dt;   % number of time steps
tcurr = 0;         % set current time step

%%  Initialize Mesh Refinement Values
% The following variables define the number of elements and their
% respective nodes in refinement 'zones.'  The total length of the sample
% is also defined.
L = 0.037;      % Case 7


% Define refinement zones with first value being 0 and last being L
%   NOTE: The refinement zone divisions MUST line up with the material
%   divisions defined below
mesh = [0 0.006 0.012 L];   % Case 7, 8, 9

% Define the number of elements in each division
%   NOTE: Length of 'ne' matrix should be one less than 'mesh' matrix
mesh_ne = [20 10 10];    % Case 7, 8, 9

% Calculate the total number of elements & nodes in the material
ne = 0;                         % initialize number of elements variable
for i = 1:length(mesh_ne)
    ne = ne + mesh_ne(i);       
end
n = ne + 1;                     % for 1D: nodes = number of elements + 1

% Calculate the node spacing through the model
mesh_dx = zeros(1,length(mesh_ne));
for i = 1:length(mesh_ne)
    mesh_dx(i) = (mesh(i+1)-mesh(i))/mesh_ne(i);
end

% Determine individual element node spacing using node spacing calculated
% in [mesh_dx]
dx = zeros(ne,1);               % initialize dx vector
dx = Decomp1D_GenerateDX( dx, mesh_ne, mesh_dx, ne );

%%  Thermal & Physical Properties of the System & Material
% The following variables define the thermal properties of the materials
% being investigated in the code.
num_mtl = 2;    % Case 7, 8, 9

% Define the divisions of material in the system
mtl_div(1) = 0.012;       % Case 7, 8, 9
mtl_div(2) = L;           % Case 7, 8, 9

% Determine node associated with material division and add to matrix
node_div = zeros(1,num_mtl,'int32');
node_div = Decomp1D_GenerateNodeDiv( node_div, num_mtl, mtl_div, mesh, mesh_ne );

%% Material Property Definitions
% Case 7, 8, 9
matl = 1;
%%% MATERIAL - E-Glass VE Laminate %%%%%%%%%%
% The following variables are used to determine properties as a function 
% of temperature for material 1:
%   Polynomial coefficients of thermal conductivity: f(T)=ak+bk*T+ck*T^2
%       kv - virgin material (kW/m-K)
%       kc - char material (kW/m-K)
akv(matl)=akv1;
bkv(matl)=bkv1;
ckv(matl)=ckv1;
akc(matl)=1.7641e-5; 
bkc(matl)=2.830e-7; 
ckc(matl)=0.0;

%   Polynomial coefficients of specific heat capacity: f(T)=ac+bc*T+cc*T^2
%       cv - virgin material (kJ/kg-K)
%       cc - char material (kJ/kg-K)
acv(matl)=acv1; 
bcv(matl)=bcv1; 
ccv(matl)=ccv1;
acc(matl)=9.703e-1; 
bcc(matl)=2.590e-4; 
ccc(matl)=0.0;

%   Polynomial coefficients of gas specific heat: f(T)=acg+bcg*T+ccg*T^2
%       cg - specific heat of gas (kJ/kg-K)
acg(matl)=-9.1151e-2; bcg(matl)=4.4007e-3; ccg(matl)=-1.7297e-6;

% The following variables define the density of the virgin & decomposed
% materials
rhov(matl)=1683.;   % virgin density (kg/m^3)
rhod(matl)=1235.;   % decomposed density (kg/m^3)

% The following materials are constants for the Arrhenius kinetics equation
% for decomposition
decPhases(matl)=1;              % number of phases in arrhenius decomposition
% Define 1st decomposition phase
currPhase=1;                    % current phase definition
Aak(matl,currPhase)=5.0e28;     % pre-exponential factor arrhenius kinetics (s^-1)
Eak(matl,currPhase)=3.62e5;     % activation energy arrhenius kinetics (kJ/kmol)
nak(matl,currPhase)=4.6;        % density fraction power arrhenius kinetics (- -)
trans(matl,currPhase)=0.0;      % transition point (F) for current phase

% Heat of Reaction
dhrec(matl,currPhase)=-500.0; % heat of reaction (decomposition) (kJ/kg)

% Set if material is inert (1 - true; 0 - false)
inert(matl) = 0;

matl = 2;
%%% MATERIAL - Superwool 607 %%%%%%%%%%
% The following variables are used to determine properties as a function 
% of temperature for material 2:
%   Polynomial coefficients of thermal conductivity: f(T)=ak+bk*T+ck*T^2
%       kv - virgin material (kW/m-K)
%       kc - char material (kW/m-K)
akv(matl)=6.696e-7; bkv(matl)=1.369e-8; ckv(matl)=2.266e-10;
akc(matl)=6.696e-7; bkc(matl)=1.369e-8; ckc(matl)=2.266e-10;

%   Polynomial coefficients of specific heat capacity: f(T)=ac+bc*T+cc*T^2
%       cv - virgin material (kJ/kg-K)
%       cc - char material (kJ/kg-K)
acv(matl)=6.542e-1; bcv(matl)=4.795e-4; ccv(matl)=-1.173e-7;
acc(matl)=6.542e-1; bcc(matl)=4.795e-4; ccc(matl)=-1.173e-7;

%   Polynomial coefficients of gas specific heat: f(T)=acg+bcg*T+ccg*T^2
%       cg - specific heat of gas (kJ/kg-K)
acg(matl)=0.0; bcg(matl)=0.0; ccg(matl)=0.0;

% The following variables define the density of the virgin & decomposed
% materials
rhov(matl)=96.;   % virgin density (kg/m^3)
rhod(matl)=96;   % decomposed density (kg/m^3)

% The following materials are constants for the Arrhenius kinetics equation
% for decomposition
decPhases(matl)=1;              % number of phases in arrhenius decomposition
% Define 1st decomposition phase
currPhase=1;                    % current phase definition
Aak(matl,currPhase)=0;     % pre-exponential factor arrhenius kinetics (s^-1)
Eak(matl,currPhase)=0;     % activation energy arrhenius kinetics (kJ/kmol)
nak(matl,currPhase)=0;        % density fraction power arrhenius kinetics (- -)
trans(matl,currPhase)=0.0;      % transition point (F) for current phase

% Heat of Reaction
dhrec(matl,currPhase)=-500.0; % heat of reaction (decomposition) (kJ/kg)

% Set if material is inert (1 - true; 0 - false)
inert(matl) = 1;
%%  Additional Properties of the System & Material

Tref=273.;    % reference temperature (K)

Rbar=8.314;   % ideal gas constant (kJ/kmol-K)
Q=0;          % heat generation (kW/m^3)

%%  Boundary Conditions

% Radiation Constants
sigma   = 5.67e-11; % Stefan-Boltzman constant (kW/m2-K4)
e1      = 0.0;      % emissivity of surface at node1 (e1=0 for no radiation)
e2      = 0.0;      % emissivity of surface at node2 (e2=0 for no radiation)

% Heat Flux & Heat Transfer Constants
%   q1/n - defined as positive into the material
%   NOTE: Each row of the following matrices is defined as follows: [time
%   q"].  The heat flux at each time step is linearly interpolated using
%   these vectors.
q1vector = [0 25; tend 25];   % Case 7
qnvector = [0 0; tend 0];     % Case 1-15

ho1      = 0.01;   % heat transfer coefficient at node1 (kW/m2-K) Case 1-15
hon      = 0.01;   % heat transfer coefficient at noden (kW/m2-K) Case 1-23

% The following variables are used to determine the external gas
% temperatures as a function of time:
%   Polynomial coefficients of To1: f(t)=aT+bT*t+cT*t^2+dT*t^3+eT*t^4
% Case 1-6
aTo1=300.0; bTo1=0.0; cTo1=0.0; dTo1=0.0; eTo1=0.0;

% Polynomial coefficients of Ton: f(t)=aT+bT*t+cT*t^2+dT*t^3+eT*t^4
% aTon=300.0; bTon=0.0; cTon=0.0; dTon=0.0; eTon=0.0;
aTon=300.0; bTon=0.0; cTon=0.0; dTon=0.0; eTon=0.0;     % Case 1-23

% Initialize external gas temperatures
To1     = eTo1*tcurr^4+dTo1*tcurr^3+cTo1*tcurr^2+bTo1*tcurr+aTo1; % Case 1-15
Ton     = eTon*tcurr^4+dTon*tcurr^3+cTon*tcurr^2+bTon*tcurr+aTon;

% Initialize internal temperatures for the sample
%  NOTE: Each row of the initial material temperature matrix is defined aso1
%  follows: [location temperature].  The initial temperature at each node
%  is lineraly interpolated using this matrrix.
% Case 7,8,9
Tinitvector = [0 300.0; 0.006 300.0; 0.012 300.0; L 300.0];

Tinit = Decomp1D_InterpolateTinit( Tinitvector, dx, n );    % initial temperature (K)
% Assign isothermal temperatures (if applicable)
T1      = 0.0;      % temperature at node1 (isothermal>0) (K)
Tn      = 0.0;      % temperature at noden (isothermal>0) (K)

% Mass Flux at Node (n)
mfn     = 0.0;      % mass flux at noden  (kg/(s-m^2) 

%%  Residuals and Iterations Control

% Define the maximum residual error and initialize the residual matrix
er      = 0.0001;        % residual error
resid   = ones(n,1)*er;  % residual matrix

% Define the maximum number of iterations and initialize the iteration
% counter
maxiter = 50;   % maximum iterations for convergence at a single time step
icount  = 0;    % iteration counter for convergence

% The following is used to define the calculation technique
theta   = 1.00;     % transient relaxation parameter

%%  Apply Initial Conditions
% This portion of code initializes all of the variables to their respective
% initial conditions for use in the "Main Time Step" loop

% Initialize the x-dimension & time step matrices
xdim      = zeros(n+(num_mtl-1),1); % x-dimension matrix for calculations
xdim_plot = ones(n,1);              % x-dimension matrix for plots
tdim      = ones(nt+1,1);           % time step matrix

% Initialize temperature values for the main time loop
Told    = ones(n,1);        % intermediate temperature in variable 
                                % properties
T       = Tinit;            % initial temperature
Tp      = T;                % temperature at previous time step
Tall    = T;                % saved temperature for time zero
Tr      = ones(ne,1)*Tref;  % reference temperature

% Calculate initial average elemental temperatures from initial temps
Tea1 = T;  Tea1(n) = [];  % set the node 1 temps; clear the nth node
Tea2 = T;  Tea2(1) = [];  % set the node 2 temps; clear the 1st node
Teavg   = (Tea1+Tea2)/2;  % initial average temperature for each element

% Define temperature profile sampling points and initialize
%Tpro_a = 1; Tpro_b = 21; Tpro_c = 31; Tpro_d = 41;    % node sampling pts for temp profile
%Tpro    = [Tinit(Tpro_a), Tinit(Tpro_b), Tinit(Tpro_c), ...
%    Tinit(Tpro_d), Tinit(Tpro_e)];  % temperature profile with time

% Dynamically builds the Tpro Matrix so that the values don't have to be
% added by hand
temp_matrix_i = 1;
temp_matrix = zeros(1,length(mesh_ne)+1);
temp_matrix(1) = Tinit(temp_matrix_i);
for t_i=1:length(mesh_ne)    
    temp_matrix_i = temp_matrix_i + mesh_ne(t_i);
    temp_matrix(t_i+1) = Tinit(temp_matrix_i);
end

Tpro = temp_matrix;

% Initialize the mass & density properties of the material
mf      = zeros(n+(num_mtl-1),1);   % intermediate mass flux
mf(n)   = mfn;                      % mass flux at unexposed side
mfp     = mf;                       % gas mass flux at previous time zero
mfpro   = zeros(1,3);               % mass flux profile with time
mfsurf  = 0;                        % mass flux at surface at time zero (kg/m2-s)

% Initialize the density properties of the material 
%   NOTE: Densities are dually defined at material intersections
drdt    = zeros(n+(num_mtl-1),1);     % change in density with time (kg/s)
r       = ones(n+(num_mtl-1),1);      % virgin material density matrix
% Assign the virgin densities to the density matrix using a subroutine
r       = Decomp1D_InitializeDensity( r, n, num_mtl, node_div, rhov );
rp      = r;                % previous virgin material density matrix
drho    = rhov-rhod;        % difference of virgin and decomposed densities

% Initialize the instantaneous mass fraction terms
Ffrac   = ones(n+(num_mtl-1),1);% instantaneous mass fraction
Feavg   = ones(ne,1);           % average mass fraction across element

% Initialize the thermal properties matrices 
kv  = zeros(ne,1);          % virgin conductivity (kW/m-K)
kc  = zeros(ne,1);          % char conductivity (kW/m-K)
cv  = zeros(ne,1);          % virgin specific heat (kJ/kg-K)
cc  = zeros(ne,1);          % char specific heat (kJ/kg-K)
cg  = zeros(ne,1);          % specific heat of gas (kJ/kg-K)
dh  = zeros(ne,1);          % enthalpy
dhg = zeros(ne,1);          % enthalpy of gas

% Calculate temperature dependent constants for all materials
kv = Decomp1D_AssignThermalProps( kv, ne, Teavg, node_div, akv, bkv, ckv );
kc = Decomp1D_AssignThermalProps( kc, ne, Teavg, node_div, akc, bkc, ckc );
cv = Decomp1D_AssignThermalProps( cv, ne, Teavg, node_div, acv, bcv, ccv );
cc = Decomp1D_AssignThermalProps( cc, ne, Teavg, node_div, acc, bcc, ccc );
cg = Decomp1D_AssignThermalProps( cg, ne, Teavg, node_div, acg, bcg, ccg );
dh = Decomp1D_AssignEnthalpies( dh, ne, Teavg, Tr, node_div, acv, bcv, ccv, ...
    dhrec, trans, Feavg );
dhg = Decomp1D_AssignEnthalpies( dhg, ne, Teavg, Tr, node_div, acg, bcg, ccg, ...
    zeros(num_mtl,num_mtl), trans, Feavg );

% Calculate the thermal properties as a function of decomposed material
k   = Feavg.*kv+(1-Feavg).*kc;        % thermal conductivity (kW/m-K)
c   = Feavg.*cv+(1-Feavg).*cc;        % specfic heat (kJ/kg-K)

% Initialize the heat transfer due to radiation
hr1 = e1*sigma*(To1+Tinit(1))*(To1^2+Tinit(1)^2); % radiation at node 1
hrn = e2*sigma*(Ton+Tinit(n))*(Ton^2+Tinit(n)^2); % radiation at node n

% Define the total heat transfer coefficient at the first and last nodes
h1  = ho1+hr1;   % total heat transfer coefficient at node 1
hn  = hon+hrn;   % total heat transfer coefficient at node n

% Generate the x-dimension matrices
xdim = Decomp1D_GenerateXDim( xdim, num_mtl, node_div, n, dx );
xdim_plot = Decomp1D_GenerateXDimPlot( xdim_plot, n, dx );

% Initialize the initial time step to zero
  tdim(1)=0; 
  
% Generate stiffness matrix
  Ap = zeros(n);    % Initialize the matrix to the number of nodes
  % This subroutine adds the heat transfer coefficients (conduction) to the
  % stiffness matrix to create a symmetric diagonal matrix:
  %     k - thermal conductivity
  %     h1 - convection coefficient at node 1
  %     hn - convection coefficient an node n
  %     dx - node spacing
  Ap = Decomp1D_AssembleStiffnessMatrix( Ap, k, h1, hn, dx, ne );

% Generate Transient Matrix (mass matrix)
  Bp = zeros(n);    % Initialize the matrix to the number of nodes
  % This loop adds the terms to the [B] matrix:
  %     r - virgin density matrix
  %     c - specific heat of the material
  %     dx - node spacing
  Bp = Decomp1D_AssembleTransientMatrix( Bp, r, c, dx, node_div, ne );

% Generate density versus time term matrix
  EDTp = zeros(n+(num_mtl-1));  % Initialize the matrix to the number of nodes
  % This loop adds the appropriate terms to the [EDT] matrix:
  %     dh - enthalpy
  %     dhg - enthalpy of gas
  %     dx - node spacing
  EDTp = Decomp1D_AssembleEDTMatrix( EDTp, dh, dhg, dx, node_div, ne );

% Generate convective term matrix
  CONVp = zeros(n); % Initialize the matrix to the number of nodes
  % This loop generates the transient heat transfer characteristics of the
  % elements for the [CONV] matrix:
  %     mf - intermediate mass flux
  %     cg - specific heat of the gas
  CONVp = Decomp1D_AssembleConvectionMatrix( CONVp, mf, cg, node_div, ne );

% Generate forcing function matrix
  Fp=zeros(n,1);    % Initialize a column vector to the number of nodes
  % This loop creates the forcing function for the equation using heat
  % generation and heat flux (for the exposed nodes)
  %     dx - node spacing
  %     Q - heat generation (separated to limit one term at nodes 1&n)
  %     q1/n - heat flux at node 1/n
  %     h1/n - heat transfer coeff at node 1/n
  %     To1/n - gas temperature at node 1/n
  [q1, qn] = Decomp1D_InterpolateHeatFlux( q1vector, qnvector, tcurr ); % Case 1-15
  Fp = Decomp1D_AssembleForcingFunction( Fp, Q, q1, qn, h1, hn, To1, Ton, dx, ne );  

%% Main Time Step Loop
% This loop iterates through the defined equation until a set number of
% times steps has been reached.  The temperatures are converged at each
% individual time step
count = 1;
dispctr = 1;
% Initialize the loop for the number of defined time steps
for j=1:nt
   
% Reinitialize and create an initial calculation for the maximum
% temperature difference at all nodes
Told     = ones(n,1);   % reinitialize the old temperature matrix from previous loop
Tdiff    = abs(Told-T); % determines the diefference between old temps and new
Tmaxdiff = max(Tdiff);  % determine the maximum temperature diff at all nodes

tcurr       = j*dt;     % set current time step to time interval * count
tdim(j+1)   = tcurr;    % set the time of the next time step as the end of 
                            % the current time step
icount      = 0;        % reinitialize the iteration ctr for the current 
                            % time step loop
% tdimstr(j+1)= {'t='}+num2cell(tcurr);

% The following loop is the "Property Update" loop, which continues to
% update the properties at each time step until the temperatures converge
% to a preset value
while Tmaxdiff > er 
  
  icount = icount+1;    % increment the iteration counter
  Told   = T;           % set the old temperatures from the previous loop
  
  % Reinitialize external gas temperatures at current time step
  To1     = eTo1*tcurr^4+dTo1*tcurr^3+cTo1*tcurr^2+bTo1*tcurr+aTo1; % Case 1-15
  Ton = eTon*tcurr^4+dTon*tcurr^3+cTon*tcurr^2+bTon*tcurr+aTon;

  % Calculate the average elemental temperature
  Tea1 = T;  Tea1(n) = [];  % set the node 1 temps; clear the nth node
  Tea2 = T;  Tea2(1) = [];  % set the node 2 temps; clear the 1st node
  Teavg = (Tea1+Tea2)/2;    % calculate the elemental avg temps
  
  % Calculate the average elemental mass fractions
  Fea1 = Ffrac;  Fea1(n) = [];  % set node 1 fractions; clear the nth node
  Fea2 = Ffrac;  Fea2(1) = [];  % set node 2 fractions; clear the 1st node
  Feavg = Decomp1D_CalculateFeavg( Fea1, Fea2, node_div, n );

  % Calculate the thermal properties for each element using the elemental
  % average temperatures.  These are lumped properties.
  %     kv/c - virgin/char thermal conductivity
  %     cv/c - virgin/char specific heat of the material
  %     cg - specific heat of the gas
  %     dh - enthalpy of the material
  %     dhg - enthalpy of the gas
  kv = Decomp1D_AssignThermalProps( kv, ne, Teavg, node_div, akv, bkv, ckv );
  kc = Decomp1D_AssignThermalProps( kc, ne, Teavg, node_div, akc, bkc, ckc );
  cv = Decomp1D_AssignThermalProps( cv, ne, Teavg, node_div, acv, bcv, ccv );
  cc = Decomp1D_AssignThermalProps( cc, ne, Teavg, node_div, acc, bcc, ccc );
  cg = Decomp1D_AssignThermalProps( cg, ne, Teavg, node_div, acg, bcg, ccg );
  dh = Decomp1D_AssignEnthalpies( dh, ne, Teavg, Tr, node_div, acv, bcv, ccv, ...
      dhrec, trans, Feavg );
  dhg = Decomp1D_AssignEnthalpies( dhg, ne, Teavg, Tr, node_div, acg, bcg, ccg, ...
      zeros(num_mtl,num_mtl), trans, Feavg );

  % Calculate the thermal properties as a function of decomposing material
  %     k - thermal conductivity
  %     c - specific heat of material
  k   = Feavg.*kv+(1-Feavg).*kc;
  c   = Feavg.*cv+(1-Feavg).*cc;

  % Calculate the heat transfer coefficients using the exterior gas
  % temperatues (To1/n) and the temperatures calculated at the previous
  % time step (Told1/n)
  hr1   = e1*sigma*(To1+Told(1))*(To1^2+Told(1)^2);  % radiation at node 1
  hrn   = e2*sigma*(Ton+Told(n))*(Ton^2+Told(n)^2);  % radiation at node n
  h1    = ho1+hr1;   % total heat transfer coeff at node 1
  hn    = hon+hrn;   % total heat transfer coeff at node n

% Generate stiffness matrix
  A=zeros(n);   % Initialize the matrix to the number of nodes
  % This loop adds the heat transfer coefficients (conduction) to the
  % stiffness matrix to create a symmetric diagonal matrix:
  %     k - thermal conductivity
  %     dx - node spacing
  A = Decomp1D_AssembleStiffnessMatrix( A, k, h1, hn, dx, ne );
% Generate Transient Matrix (mass matrix)
  B=zeros(n);   % Initialize the matrix to the number of nodes
  % This loop adds the terms to the [B] matrix:
  %     r - virgin density matrix
  %     c - specific heat of the material
  %     dx - node spacing
  B = Decomp1D_AssembleTransientMatrix( B, r, c, dx, node_div, ne );

% Generate density versus time term matrix
  EDT=zeros(n+(num_mtl-1)); % Initialize the matrix to the number of nodes
  % This loop adds the appropriate terms to the [EDT] matrix:
  %     dh - enthalpy
  %     dhg - enthalpy of gas
  %     dx - node spacing
  EDT = Decomp1D_AssembleEDTMatrix( EDT, dh, dhg, dx, node_div, ne );
 
 
% Generate convective term matrix
  CONV=zeros(n);    % Initialize the matrix to the number of nodes
  % This loop generates the transient heat transfer characteristics of the
  % elements for the [CONV] matrix:
  %     mf - intermediate mass flux
  %     cg - specific heat of the gas
  CONV = Decomp1D_AssembleConvectionMatrix( CONV, mf, cg, node_div, ne );   

  %Generate forcing function matrix
  F=zeros(n,1);     % Initialize a column vector to the number of nodes
  % This loop creates the forcing function for the equation using heat
  % generation and heat flux (for the exposed nodes)
  %     dx - node spacing
  %     Q - heat generation (separated to limit one term at nodes 1&n)
  %     q1/n - heat flux at node 1/n
  %     h1/n - heat transfer coeff at node 1/n
  %     To1/n - gas temperature at node 1/n
  [q1, qn] = Decomp1D_InterpolateHeatFlux( q1vector, qnvector, tcurr ); % Case 1-15
  F = Decomp1D_AssembleForcingFunction( F, Q, q1, qn, h1, hn, To1, Ton, dx, ne );   

  %Premultiply the EDT and drdt matrices
  EDT_drdt=zeros(n,1);  % Initiliaze the matrix to the number of nodes
  % This loop premultiplies the EDT and drdt matrices in a manner which
  % allows the use of multiple materials
  EDT_drdt = Decomp1D_PremultiplyEDTdrdt( EDT_drdt, EDT, drdt, node_div, n );
  
% GENERATE MAIN MATRICES

  % Matrix with current time step properties
  D = B+theta*dt*(A+CONV);
  % Matrix with previous time step values
  E = B-(1-theta)*dt*(A+CONV);
  % Forcing matrix for current time step
  G = E*Tp+(1-theta)*dt*Fp+theta*dt*F-EDT_drdt*dt;

% Isothermal BC - Set temperature to isothermal temperature (Constant
% Temperature BC)
  if (T1>0) 
     D(1,:) = zeros;
     D(1,1) = 1;
     G(1)   = T1;
  end 

  if (Tn>0) 
     D(n,:) = zeros;
     D(n,n) = 1;
     G(n)   = Tn;
  end 

% Calculate the Temperature Distribution through all Elements
  T = D\G;    % Determined through manipulation of orginial equation
   
  Tdiff     = abs(Told-T);  % previous time step temperature difference
  Tmaxdiff  = max(Tdiff);   % max temperature difference
 
  % Check to is if maximum iterations have been reach

end % End the "Property Update" loop

%%%%%%%%%%%%%%%%%%
% Density and Density change with Time
% Mass Conservation
%%%%%%%%%%%%%%%%%%

  % Calculate the density change with time for all elements
  % Terms:
  %     r - material density matrix
  %     drdt - density change with time
  %     drho - difference of virgin and decomposed densities
  %     rhod - decomposed density
  %     nak - density fraction power - Arrhenius kinetics
  %     Aak - pre-exponential factor - Arrhenius kinetics
  %     Eak - activation energy - Arrhenius kinetics
  %     Rbar - universal gas constant
  %     T(i) - temperature for node i
  drdt = Decomp1D_CalculateDRDT( drdt, r, drho, rhod, nak, Aak, Eak,...
      trans, Ffrac, Rbar, T, inert, n, node_div );
  
  % Mass flux calculation
  mf = Decomp1D_CalculateMassFlux( mf, xdim, drdt, n, ne, num_mtl );
  
  % Density updates - r(new) = dr/dt * dt + r(old)
  rp = r;           % set previous virgin density to be current density
  r  = drdt*dt+rp;  % calculate the new density matrix

  % Instantaneous mass fraction updates
  Ffrac = Decomp1D_UpdateMassFraction( Ffrac, r, rhod, drho, n, node_div,count );
  count = count+1;
% Horizontally concatenate temperatures in one matrix
Tall    = horzcat(Tall,T);
% Vertically concatenate mass fluxes into one matrix
mfsurf  = vertcat(mfsurf,-1*mf(1));

% Dynamically builds the Tpro Matrix so that the values don't have to be
% added by hand
temp_matrix_i = 1;
temp_matrix = zeros(1,length(mesh_ne)+1);
temp_matrix(1) = T(temp_matrix_i);
for t_i=1:length(mesh_ne)
    temp_matrix_i = temp_matrix_i + mesh_ne(t_i);
    temp_matrix(t_i+1) = T(temp_matrix_i);
end

% Vertically concatenate the temperatures at the first, middle, and last
% nodes into a single matrix
Tpro    = vertcat(Tpro,temp_matrix);

% Assign current solution to be solution at previous time step for the next
% iteration
Ap      = A;
Bp      = B;
EDTp    = EDT;
CONVp   = CONV;
Fp      = F;
Tp      = T;
mfp     = mf;

end  % End the "Main Time Step" loop