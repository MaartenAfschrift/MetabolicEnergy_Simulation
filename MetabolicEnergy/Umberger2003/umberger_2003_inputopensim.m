function [ musclemass, slowtwitchratio, maximalcontractionvelocity] = umberger_2003_inputopensim(modelpath, muscle)
% UMBERGERINPUTOSIM
% parameters needed to calculate metabolic energy consumption for a muscle
% are derived from opensim

import org.opensim.modeling.*
mymodel = Model(modelpath);
mystate = mymodel.initSystem();
mymuscleset = mymodel.getMuscles();
mymuscle = mymuscleset.get(muscle);
mymuscleprobe = Umberger2010MuscleMetabolicsProbe();
mymuscleprobe.addMuscle(muscle, 0.70);
mymodel.addProbe(mymuscleprobe);

musclemass = mymuscleprobe.getMuscleMass(muscle);
% tension = mymuscleprobe.getSpecificTension(muscle);
slowtwitchratio = mymuscleprobe.getRatioSlowTwitchFibers(muscle);
% density = mymuscleprobe.getDensity(muscle);
% optimalfiberlength = mymuscle.get_optimal_fiber_length();
% optimalfiberforce = mymuscle.get_max_isometric_force();
maximalcontractionvelocity = mymuscle.get_max_contraction_velocity();
% control_musclemass = density*optimalfiberforce*optimalfiberlength/tension;

end