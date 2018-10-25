% This function computes the muscle energy expenditure based on the model
% of Umberger2010. Energy_total is in W. 
% This function provides a smooth approximation of the original function
% for use in optimization problems.
% As compared to getMetabolicEnergySmoothUmberger2010, this versions 1)
% includes the basal heat rate and 2) insures that the total heat rate is
% at least equal to 1.
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
%   where rho is density (1058.7 kg/m�) and sigma is tension (0.25MPa)
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

function [hdottotal, energy_am, energy_sl] = calculateheat_2010(exc,act,lMtilde,vMtilde,pctst,vcemax,Fiso,b)

%% Some parameters
% Percentage of slow twitch (st) and fast twitch (ft) fibers
pctft = 1-pctst;
% Scaling factor: aerobic (1.5) or anaerobic (1) activities
s = 1.5;            

%% Excitations and activations
a = exc + (-exc+act)/2.*(0.5+0.5*tanh(b*(act-exc)));

%% Heat activation and maintenanance
% lMtilde is on the descending leg of the active force-length relationship
% i.e. lMtilde > 1 || lM > lMopt
lMtilde_des = 0.5 + 0.5*tanh(b*(lMtilde-1*ones(size(lMtilde,1),1)));        
hdotam = 128*pctft+25 + (-0.6*(128*pctft+25) + 0.6*(128*pctft+25) ... 
         .*Fiso).*lMtilde_des;       
aam         = a.^(0.6);
energy_am   = hdotam.*aam*s;

%% Heat shortening and lengthening
% vM is positive (muscle lengthening)
vMtilde_pos = 0.5 + 0.5*tanh(b*(vMtilde));
% vM is negative (muscle shortening)
vMtilde_neg = 1-vMtilde_pos;

coef_hs_ft  = 1*153.*ones(size(vcemax,1),1)./vcemax;
coef_hs_st  = 4*25.*ones(size(vcemax,1),1)./(vcemax/2.5); 
coef_hl     = 0.3*coef_hs_st; % different as compared to Umberger2003

hdotsl = coef_hl.*vMtilde + (-coef_hl.*vMtilde ...
    - coef_hs_st.*vMtilde.*pctst ...
    - coef_hs_ft.*vMtilde.*pctft).*vMtilde_neg;
energy_sl = (hdotsl.*a.*s + ...
    (-hdotsl.*a.*s + hdotsl.*(a.^2)*s).*vMtilde_neg).*(1-lMtilde_des) + ...
    (hdotsl.*a.*s.*Fiso + (-hdotsl.*a.*s.*Fiso + ...
    hdotsl.*(a.^2)*s.*Fiso).*vMtilde_neg).*(lMtilde_des);

%% Heat total 

hdottotal  = energy_am + energy_sl;
end

