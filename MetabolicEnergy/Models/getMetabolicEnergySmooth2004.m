% This function computes the muscle energy expenditure based on the model
% of Bhargava2004.
%
% TODO: no fiber length dependence on maintenance rate for now. We only
% use the fiber length. In OpenSim, they use a Piecewise Linear Function
% https://github.com/opensim-org/opensim-core/blob/master/OpenSim/Simulation/Model/Bhargava2004MuscleMetabolicsProbe.cpp#L103
%
% Note: we consider some lenghtening heat rate, i.e. 
% get_use_force_dependent_shortening_prop_constant is true (OpenSim)
%
% Date: 25/10/2018
% Inputs:
%   exc: muscle excitation
%   act: muscle activation
%   lMtilde: normalized muscle fiber length
%   vM: muscle fiber velocity (positive is lengthening)
%   Fce: muscle force from contractile element (length + velocity
%   components but not passive component)
%   Fpass: passive muscle focre
%   musclemass: mass of muscle, mass = PCSA*rho*lMopt, PCSA = Fmax/sigma
%   where rho is density (1058.7 kg/m³) and sigma is tension (0.25MPa)
%   pctst: percentage of slow twitch fibers (0-1)
%   vcemax: maximal muscle fiber velocity (default is 10*lMopt)
%   Fiso: normalized muscle force from active force-length relationship
% Outputs:
%   energy_total: total energy rate
%   Adot: energy rate from activation
%   Mdot: energy rate from maintenance
%   Sdot: energy rate from shortening and lengthening
%   Wdot: energy rate from mechanical work 
%   energy_model: energy rate model including basal rate
% Note: the energy_total might be different than the sum of the different
% component if the total heat rate was clamped to one.

function [energy_total,Adot,Mdot,Sdot,Wdot,energy_model] = ...
    getMetabolicEnergySmooth2004(exc,act,lMtilde,vM,Fce,Fpass,...
                                musclemass,pctst,Fiso,Fmax,modelmass,b)

%% Some parameters
% Percentage of slow twitch (st) and fast twitch (ft) fibers
pctft = 1-pctst;      

%% Twitch excitations
st_e = pctst.*sin(pi/2*exc);
ft_e = pctft.*(1-cos(pi/2*exc));

%% Activation heat rate
decay_function_value = 1;
activation_constant_st = 40; % default (Bhargava et al. 2004)
activation_constant_ft = 133; % default (Bhargava et al. 2004)
Adot = musclemass.*decay_function_value.*((activation_constant_st*st_e)+...
    (activation_constant_ft*ft_e));

%% Maintenace heat rate
fiber_length_dep = lMtilde; % TODO OpenSim uses a Piecewise Linear Function
maintenance_constant_st = 74; % default (Bhargava et al. 2004)
maintenance_constant_ft = 111; % default (Bhargava et al. 2004)
Mdot = musclemass.*fiber_length_dep.*((maintenance_constant_st*st_e)+...
    (maintenance_constant_ft*ft_e));

%% Shortening heat rate
% F_iso that would be developed at the current activation and fiber length
% under isometric conditions (THIS DIFFERS FROM UMBERGER/UCHIDA)
% Fiso is getActiveForceLengthMultiplier in OpenSim. To minimize the
% difference between the models, we keep the same input, i.e. Fiso, that is
% used in Umberger/Uchida and here.
F_iso = act.*Fiso.*Fmax; 
fiber_force_total = Fce + Fpass;
alpha = (0.16 * F_iso) + (0.18 * fiber_force_total);
% vM is positive (muscle lengthening)
vM_pos = 0.5 + 0.5*tanh(b*(vM));
% vM is negative (muscle shortening)
vM_neg = 1-vM_pos;
alpha = alpha + (-alpha + 0.157*fiber_force_total).*vM_pos;
% If get_use_force_dependent_shortening_prop_constant is false then
% alpha = 0.25*fiber_force_total;
% alpha(vM>0) = 0;
Sdot = -alpha.*vM;

%% Mechanical work rate
% No negative mechanical work allowed
Wdot = - Fce.*vM.*vM_neg;

%% Total power is non-negative
% If necessary, increase the shortnening heat rate
% https://github.com/opensim-org/opensim-core/blob/master/OpenSim/Simulation/Model/Bhargava2004MuscleMetabolicsProbe.cpp#L393
Edot_W_beforeClamp = Adot + Mdot + Sdot + Wdot;
Edot_Wkg_beforeClamp_neg = 0.5+(0.5*tanh(b*(-Edot_W_beforeClamp)));
Sdot = Sdot - Edot_W_beforeClamp.*Edot_Wkg_beforeClamp_neg;

%% The total heat rate is not allowed to drop below 1(W/kg)
% https://github.com/opensim-org/opensim-core/blob/master/OpenSim/Simulation/Model/Bhargava2004MuscleMetabolicsProbe.cpp#L400
totalHeatRate = Adot + Mdot + Sdot;
% We first express in W/kg 
totalHeatRate = totalHeatRate./musclemass;
totalHeatRate = totalHeatRate + (-totalHeatRate + ...
    ones(size(totalHeatRate,1),1)).*(0.5 + 0.5*tanh(b*(1-totalHeatRate)));
% We then express back in W
totalHeatRate = totalHeatRate.*musclemass;

%% Total metabolic energy rate
energy_total = totalHeatRate + Wdot;

%% Energy model
basal_coef = 1.2; % default in OpenSim
basal_exp = 1; % default in OpenSim
energy_model = basal_coef*modelmass^basal_exp + sum(energy_total); 
end
