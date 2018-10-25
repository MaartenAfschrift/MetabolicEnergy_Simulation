clear all; close all; clc

%% Example_MuscleRedundacySolver
% This script is used to solve the muscle redundancy problem in a simple
% example and to simulate muscle-tendon dynamics

%% Requirements

% - Installation GPOPS II for muscle redundancy solver
% - Installation Adigator
% - OpenSim 3.2 or 3.3 for muscle analysis (including matlab API)


%% path

addpath(genpath('C:\Users\u0098084\Documents\MATLAB\Metabolic Energy'));

%% Setting

S.GetOptInput=1;
S.OpenSimInstallation='C:\OpenSim 3.3';
S.OutFolder='C:\tempData';              % please don't add temporary simulations results to the main project folder
S.BoolLinearSpring = 1;                 % using linear or non-linear spring for tendon

%% Example metabolic energy

% To Do: Add different continuous functions for linear and non-linear
% tendons. Doesn't work with logical operators due to Adigator

if S.GetOptInput==1    
    % get walking example OpenSim results
    Datapath=fullfile(S.OpenSimInstallation, 'Models\Gait10dof18musc\OutputReference');
    IK_path=fullfile(Datapath,'IK','subject01_walk_IK.mot');
    ID_path=[];%fullfile(Datapath,'ID','inversedynamics.sto'); % compute ID from the external loads
    model_path=fullfile(Datapath,'subject01.osim');
    time=[0.6 1.7];     % Part of the right stance phase
    if ~isdir(S.OutFolder ); mkdir(S.OutFolder ); end
    Misc.DofNames_Input={'ankle_angle_r'};
    Misc.Loads_path=fullfile(Datapath,'ExperimentalData','subject01_walk_grf.xml');
    Misc.ID_ResultsPath=fullfile(Datapath,'ID','inversedynamics.sto');
    Misc.MuscleNames_Input={'soleus_r'};      % Selects all muscles for the Input DOFS when this is left empty.
    Misc.Atendon=35;
    Misc.Mesh_Frequency=150;         
    [setup,auxdata,DatStore]=GetFtildeSetup(model_path,IK_path,ID_path,time,S.OutFolder ,Misc);  
    save('optInput.mat','setup','auxdata','DatStore');
else
    load('optInput.mat')
end


%% Run optimization for default tendon stiffness
ATendon = 35;
setup.auxdata.Atendon=ATendon;    
output = gpops2(setup);
res=output.result.solution.phase(1);
Time=res.time;
act=res.state(:,1:auxdata.NMuscles);
TForcetilde=res.state(:,auxdata.NMuscles+1:auxdata.NMuscles*2);
fce=TForcetilde.*(ones(size(Time))*DatStore.Fiso);
exc=res.control(:,1:auxdata.NMuscles);
RActivation=res.control(:,auxdata.NMuscles+1:auxdata.NMuscles+auxdata.Ndof);
TForce_dottilde=res.control(:,auxdata.NMuscles+auxdata.Ndof+1:end);
MuscleNames=DatStore.MuscleNames;
OptInfo=output;

lMTinterp=zeros(length(Time),auxdata.NMuscles);
VMTinterp=zeros(length(Time),auxdata.NMuscles);
for m=1:auxdata.NMuscles
    [lMTinterp(:,m),VMTinterp(:,m),~] = SplineEval_ppuval(auxdata.LMTSpline(m),Time,1);
end
[lM,lMtilde,vM,vMtilde ] = FiberVelocity_Ftilde(TForcetilde,TForce_dottilde,auxdata.params,lMTinterp,VMTinterp,ATendon,S.BoolLinearSpring);
Fiso = get_Flm_tilde(lMtilde,auxdata.Faparam);


%% compute metabolic energy consumption using different models
b = 1000;
modelmass = 75.1646;
musclemass = (5137/250000).*0.05.*1059.7; % estimation for soleus based on OpenSim and Umberger 2003 values: (optimalforce(OpenSim)/specifictension(Umberger))*optimalfiberlength(OpenSim)*muscledensity(Umberger)
pctst = 0.70; % data from Edgerton 1975 for the soleus
vcemax = 10; % max_contraction_velocity from OpenSim soleus
[energy_total03,energy_am03,energy_sl03,energy_mech03,energy_model03] = getMetabolicEnergySmooth2003all(exc,act,lMtilde,vMtilde,vM,fce,musclemass,pctst,vcemax,Fiso,modelmass,b);
[energy_total10,energy_am10,energy_sl10,energy_mech10,energy_model10] = getMetabolicEnergySmooth2010all(exc,act,lMtilde,vMtilde,vM,fce,musclemass,pctst,vcemax,Fiso,modelmass,b);
[energy_total16,energy_am16,energy_sl16,energy_mech16,energy_model16] = getMetabolicEnergySmooth2016all(exc,act,lMtilde,vMtilde,vM,fce,musclemass,pctst,vcemax,Fiso,modelmass,b);

%% plot

figure()
plot(energy_total03(5:end-5,:)); hold on; plot(energy_total10(5:end-5,:)); hold on; plot(energy_total16(5:end-5,:));
legend('Umberger 2003', 'Umberger 2010', 'Uchida 2016');
title('total energy cost');

figure()
plot(energy_am03(5:end-5,:)); hold on; plot(energy_am10(5:end-5,:)); hold on; plot(energy_am16(5:end-5,:));
legend('Umberger 2003', 'Umberger 2010', 'Uchida 2016');
title('activation and maintenance energy cost ');

figure()
plot(energy_sl03(5:end-5,:)); hold on; plot(energy_sl10(5:end-5,:)); hold on; plot(energy_sl16(5:end-5,:));
legend('Umberger 2003', 'Umberger 2010', 'Uchida 2016');
title('fiber shortening and lengthening energy cost');

figure()
plot(energy_mech03(5:end-5,:)); hold on; plot(energy_mech10(5:end-5,:)); hold on; plot(energy_mech16(5:end-5,:));
legend('Umberger 2003', 'Umberger 2010', 'Uchida 2016');
title('mechanical energy cost');
