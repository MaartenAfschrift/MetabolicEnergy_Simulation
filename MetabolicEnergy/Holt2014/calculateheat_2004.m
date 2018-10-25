% This function computes the muscle energy expenditure based on the model
% of Bhargava2004. Energy_total is in W. 
% This function provides a smooth approximation of the original function
% for use in optimization problems.
% As compared to getMetabolicEnergySmoothUmberger2010all, this versions 1)
% includes negative mechanical work, 2) ensures that the total power is 
% non-negative, 3) use an order recruitment model
%
% Code inspired from Maarten Afschrift and Azin Zargham
% Date: 14/02/2018
% Inputs:
%   act: muscle activation
%   lMtilde: normalized muscle fiber length
%   vMtilde: normalized muscle fiber velocity (positive is lengthening).
%   Note: defined for this model as vM/lMopt, which is a different metric
%   to the typical vM/vMax
%   vM: muscle fiber velocity (positive is lengthening)
%   Fce: muscle force from contractile element (length + velocity
%   components but not passive component)
%   musclemass: mass of muscle, mass = PCSA*rho*lMopt, PCSA = Fmax/sigma
%   where rho is density (1058.7 kg/m³) and sigma is tension (0.25MPa)
%   pctst: percentage of slow twitch fibers (0-1)
%   vcemax: maximal muscle fiber velocity (default is 10*lMopt)
%   Fiso: normalized muscle force from active force-length relationship
%   modelmass: mass model
%   b: parameter determining transition smoothness for tanh approximations
% Outputs:
%   energy_total: totat metabolic energy (W)
%   energy_am: energy from activation and maintenance (W)
%   energy_sl: energy from shortening and lengthening (W)
%   energy_mech: mechanical energy (W)
%   energy_model: energy model (includes basal rate) (W)
% Note: the energy_total might be different than the sum of the different
% component if the total heat rate was clamped to one.

function [hdottotal, energy_am, energy_sl] = ...
    calculateheat_2004(exc,lMtilde,vMtilde,pctst,fm)
%% Some parameters

% activation constant for ST fibers
Adot_fast=133;
% activation constant for FT fibers
Adot_slow=40;
% maintenance constant for ST fibers
Mdot_fast=111;
% maintenance constant for FT fibers
Mdot_slow=74;

%% Heat activation 

Adot = Adot_slow * pctst * sin((pi/2)* exc) + Adot_fast * (1-pctst) * (1-cos((pi/2)*exc));

%% Heat maintenanance

% normalized fiber length dependence of the maintenance heat rate
[func] = Fitted_Mdot(lMtilde);

Mdot = func .* (Mdot_slow * pctst * sin((pi/2)*exc) + Mdot_fast * (1-pctst) * (1-cos((pi/2)*exc)));

%% Heat shortening and lengthening

vp=0.5+0.5.*tanh(100.*vMtilde);
vn=1-vp;

alpha=0.25.*(fm).*(vn); % with respect to OpenSim probe alpha is calculated independent from fiber force (tick box in OpenSim probe)

Sdot= -alpha.*vMtilde;

%% Heat total

hdottotal = Adot + Mdot + Sdot;
energy_am = Adot + Mdot;
energy_sl = Sdot;

end

