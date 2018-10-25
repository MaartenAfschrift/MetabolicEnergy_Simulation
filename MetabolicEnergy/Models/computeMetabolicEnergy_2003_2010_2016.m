function [ energy ] = computeMetabolicEnergy_2003_2010_2016( input,do )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

b = 1000;
muscles = do.MuscleNames.stride1;
exc = do.MExcitation.stride1;
act = do.MActivation.stride1;
lm = do.lM.stride1;
lmtilde = do.lMtilde.stride1;
time = do.Time.stride1;
model = input.model;
import org.opensim.modeling.*;
curr_mymodel = Model(model);
curr_mystate = curr_mymodel.initSystem();
modelmass = curr_mymodel.getTotalMass(curr_mystate);
curr_muscleset = curr_mymodel.getMuscles();
for i=1:length(muscles)
    pp = spline(time,lm(:,i));
    dpp = fnder(pp,1);
    vm(:,i) = ppval(dpp,time);
    curr_muscle = curr_muscleset.get(muscles(i));
    lmopt(i) = curr_muscle.getOptimalFiberLength();
    vmmax(i) = curr_muscle.getMaxContractionVelocity();
    fmopt(i) = curr_muscle.getMaxIsometricForce();
    tension(i) = getSpecificTensions(muscles(i));
    pctst(i) = getSlowTwitchRatios(muscles(i));
    vcemaxtilde(i) = curr_muscle.getMaxContractionVelocity();
end
vcemax = vcemaxtilde.*lmopt;
lmopt_mat = ones(length(exc),1)*lmopt;
vmmax_mat = ones(length(exc),1)*vmmax;
fmopt_mat = ones(length(exc),1)*fmopt;
vmtilde_lmopt = vm./lmopt_mat; % necessary as input to compute metabolic energy consumption
vmtilde_vmmax = vm./vmmax_mat; % necessary as input to compute force-velocity
[fce,fiso] = get_fce_fiso(act,lmtilde,vmtilde_vmmax,fmopt_mat);
density = 1059.7;
musclemass = density*fmopt/tension*lmopt/1e6;

[energy.total_2003,energy.heat_am_2003,energy.heat_sl_2003,energy.mech_work_2003,energy.model_2003] = getMetabolicEnergySmooth2003allmuscles(exc,act,lmtilde,vmtilde_lmopt,vm,fce,musclemass,pctst,vcemax,fiso,b,modelmass);
[energy.total_2010,energy.heat_am_2010,energy.heat_sl_2010,energy.mech_work_2010,energy.model_2010] = getMetabolicEnergySmooth2010allmuscles(exc,act,lmtilde,vmtilde_lmopt,vm,fce,musclemass,pctst,vcemax,fiso,b,modelmass);
[energy.total_2016,energy.heat_am_2016,energy.heat_sl_2016,energy.mech_work_2016,energy.model_2016] = getMetabolicEnergySmooth2016allmuscles(exc,act,lmtilde,vmtilde_lmopt,vm,fce,musclemass,pctst,vcemax,fiso,b,modelmass);
end

