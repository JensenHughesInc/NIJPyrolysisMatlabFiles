%%Material Properties
%% MATERIAL #1
% MATERIAL - ***
matl = 1;
% The following variables are used to determine properties as a function 
% of temperature for material 1:
%   Polynomial coefficients of thermal conductivity: f(T)=ak+bk*T+ck*T^2
%       kv - virgin material (kW/m-K)
%       kc - char material (kW/m-K)
akv(matl)=2.905e-4; bkv(matl)=9.6566e-8; ckv(matl)=-5.9065e-11;
akc(matl)=1.7641e-5; bkc(matl)=2.830e-7; ckc(matl)=0.0;

%   Polynomial coefficients of specific heat capacity: f(T)=ac+bc*T+cc*T^2
%       cv - virgin material (kJ/kg-K)
%       cc - char material (kJ/kg-K)
acv(matl)=1.03; bcv(matl)=2.3061e-4; ccv(matl)=-2.0892e-7;
acc(matl)=9.703e-1; bcc(matl)=2.590e-4; ccc(matl)=0.0;

%   Polynomial coefficients of gas specific heat: f(T)=acg+bcg*T+ccg*T^2
%       cg - specific heat of gas (kJ/kg-K)
acg(matl)=-9.1151e-2; bcg(matl)=4.4007e-3; ccg(matl)=-1.7297e-6;

% The following variables define the density of the virgin & decomposed
% materials
rhov(matl)=1683.;   % virgin density (kg/m^3)
rhod(matl)=1235.;   % decomposed density (kg/m^3)
if(strcmp(SA_flag,'rhov1'))
    rhov(matl) = SA*rhov(matl);
elseif(strcmp(SA_flag,'rhod1'))
    rhod(matl) = SA*rhod(matl);
end

% The following materials are constants for the Arrhenius kinetics equation
% for decomposition
decPhases(matl)=1;              % number of phases in arrhenius decomposition
% Define 1st decomposition phase
currPhase=1;                    % current phase definition
Aak(matl,currPhase)=5.0e28;     % pre-exponential factor arrhenius kinetics (s^-1)
Eak(matl,currPhase)=3.62e5;     % activation energy arrhenius kinetics (kJ/kmol)
nak(matl,currPhase)=4.6;        % density fraction power arrhenius kinetics (- -)
trans(matl,currPhase)=0.0;      % transition point (F) for current phase

if(strcmp(SA_flag,'A'))
    Aak(matl,currPhase)=SA*Aak(matl,currPhase);
elseif(strcmp(SA_flag,'E'))
    Eak(matl,currPhase)=SA*Eak(matl,currPhase);
elseif(strcmp(SA_flag,'n'))
    nak(matl,currPhase)=SA*nak(matl,currPhase);
end

% Heat of Reaction
dhrec(matl,currPhase)=-500.0; % heat of reaction (decomposition) (kJ/kg)
if(strcmp(SA_flag,'dhrec'))
    dhrec(matl,currPhase)=SA*dhrec(matl,currPhase);
end

% Set if material is inert (1 - true; 0 - false)
inert(matl) = 0;

%% MATERIAL #2
% MATERIAL - ***
matl = 2;
% The following variables are used to determine properties as a function 
% of temperature for material 1:
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
rhod(matl)=96.;   % decomposed density (kg/m^3)

% The following materials are constants for the Arrhenius kinetics equation
% for decomposition
decPhases(matl)=1;              % number of phases in arrhenius decomposition
% Define 1st decomposition phase
currPhase=1;                    % current phase definition
Aak(matl,currPhase)=0.0;     % pre-exponential factor arrhenius kinetics (s^-1)
Eak(matl,currPhase)=0.0;     % activation energy arrhenius kinetics (kJ/kmol)
nak(matl,currPhase)=0.0;        % density fraction power arrhenius kinetics (- -)
trans(matl,currPhase)=0.0;      % transition point (F) for current phase

% Heat of Reaction
dhrec(matl,currPhase)=-500.0; % heat of reaction (decomposition) (kJ/kg)

% Set if material is inert (1 - true; 0 - false)
inert(matl) = 1;