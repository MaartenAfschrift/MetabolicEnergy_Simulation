% This function computes the muscle energy expenditure based on the model
% of Umberger2003. Energy_total is in W/kg and should be multiplied by the
% muscle mass to have the actual power. 
%
% Code from Maarten Afschrift
% Date: 14/02/2018
% Inputs:
%   exc: muscle excitation
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
% Outputs:
%   energy_total: totat metabolic energy (W/kg)
%   energy_am: energy from activation and maintenance
%   energy_sl: energy from shortening and lengthening
%   energy_mech: mechanical energy
% Note: the energy_total might be different than the sum of the different
% component if the total heat rate was clamped to one.

function [energy_total,energy_am,energy_sl,energy_mech] = ...
    getMetabolicEnergyNonSmooth2003(exc,act,lMtilde,vMtilde,vM,Fce, ...
                                musclemass,pctst,vcemax,Fiso)

%% Some parameters
% Percentage of slow twitch (st) and fast twitch (ft) fibers
pctft = 1-pctst;
% Scaling factor: aerobic (1.5) or anaerobic (1) activities
s = 1.5;            

%% Excitations and activations
ind_act_1 = find(exc > act);
ind_act_2 = find(exc <= act);

a(ind_act_1,1) = exc(ind_act_1,1);
a(ind_act_2,1) = (exc(ind_act_2,1)+act(ind_act_2,1))/2;

%% Heat activation and maintenanance
ind_lv_1 = find(lMtilde <= 1 & vMtilde <= 0);
ind_lv_2 = find(lMtilde <= 1 & vMtilde >  0);
ind_lv_3 = find(lMtilde >  1 & vMtilde <= 0);
ind_lv_4 = find(lMtilde >  1 & vMtilde >  0);

hdotam(ind_lv_1,1) = 128*pctft(ind_lv_1,1)+25;
hdotam(ind_lv_2,1) = 128*pctft(ind_lv_2,1)+25;
hdotam(ind_lv_3,1) = 0.4*(128*pctft(ind_lv_3,1)+25)+ ...
                     0.6.*(128*pctft(ind_lv_3,1)+25).*Fiso(ind_lv_3,1);
hdotam(ind_lv_4,1) = 0.4*(128*pctft(ind_lv_4,1)+25)+ ...
                     0.6.*(128*pctft(ind_lv_4,1)+25).*Fiso(ind_lv_4,1);
                 
aam = a.^(0.6);
energy_am = hdotam.*aam*s;

%% Heat shortening and lengthening
coef_hs_ft  = 1*153.*ones(size(vcemax,1),1)./vcemax;
coef_hs_st  = 4*25.*ones(size(vcemax,1),1)./(vcemax/2.5); 
coef_hl     = 4*coef_hs_st;

hdotsl(ind_lv_1,1) = -1.*coef_hs_st(ind_lv_1,1).*vMtilde(ind_lv_1,1).*pctst(ind_lv_1,1) ...
                     - coef_hs_ft(ind_lv_1,1).*vMtilde(ind_lv_1,1).*pctft(ind_lv_1,1);
hdotsl(ind_lv_2,1) = coef_hl(ind_lv_2,1).*vMtilde(ind_lv_2,1);
hdotsl(ind_lv_3,1) = -1.*coef_hs_st(ind_lv_3,1).*vMtilde(ind_lv_3,1).*pctst(ind_lv_3,1) ...
                     - coef_hs_ft(ind_lv_3,1).*vMtilde(ind_lv_3,1).*pctft(ind_lv_3,1);
hdotsl(ind_lv_4,1) = coef_hl(ind_lv_4,1).*vMtilde(ind_lv_4,1);

energy_sl(ind_lv_1,1) = hdotsl(ind_lv_1,1).*(a(ind_lv_1,1).^2).*s;
energy_sl(ind_lv_2,1) = hdotsl(ind_lv_2,1).*a(ind_lv_2,1).*s;
energy_sl(ind_lv_3,1) = hdotsl(ind_lv_3,1).*Fiso(ind_lv_3,1) ...
                        .*(a(ind_lv_3,1).^2).*s; 
energy_sl(ind_lv_4,1) = hdotsl(ind_lv_4,1).*Fiso(ind_lv_4,1) ...
                        .*a(ind_lv_4,1).*s;

%% Mechanical work
% This takes into account negative mechanical work
energy_mech = -1*Fce.*vM./musclemass;

%% The total rate of energy liberation is not allowed to drop below 1(W/kg)
% https://github.com/opensim-org/opensim-core/blob/master/OpenSim/Simulation/Model/Umberger2010MuscleMetabolicsProbe.cpp#L432
totalHeatRate = energy_am + energy_sl;
ind_thr_1 = totalHeatRate < 1;
totalHeatRate(ind_thr_1,1) = 1;  

%% Account for muscle mass
energy_am = energy_am.*musclemass;
energy_sl = energy_sl.*musclemass;
energy_mech = energy_mech.*musclemass;
energy_total = totalHeatRate.*musclemass+energy_mech;

end

