clear all; close all; clc

%% Example_MuscleRedundacySolver
% This script is used to solve the muscle redundancy problem in a simple
% example and to simulate muscle-tendon dynamics

%% Requirements

% - Installation GPOPS II for muscle redundancy solver
% - Installation Adigator
% - OpenSim 3.2 or 3.3 for muscle analysis (including matlab API)

%% Notes

% - It would be interesting to see what happens when you add a tibial
% muscle  (Maarten
% - try this example on the hamner running data as well  (Maarten)

%% Add paths

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
output = gpops2(setup);
res=output.result.solution.phase(1);
Time=res.time;
MActivation=res.state(:,1:auxdata.NMuscles);
TForcetilde=res.state(:,auxdata.NMuscles+1:auxdata.NMuscles*2);
TForce=TForcetilde.*(ones(size(Time))*DatStore.Fiso);
MExcitation=res.control(:,1:auxdata.NMuscles);
RActivation=res.control(:,auxdata.NMuscles+1:auxdata.NMuscles+auxdata.Ndof);
MuscleNames=DatStore.MuscleNames;
OptInfo=output;
lMTinterp = interp1(DatStore.time,DatStore.LMT,Time);
[lM,lMtilde] = FiberLength_Ftilde(TForcetilde,auxdata.params,lMTinterp);

%% Run optimization for different levels of tendon stiffness

% Compute energy-related information on the skeletal level
Ta=ppval(auxdata.JointIDSpline,Time);
IK=importdata(IK_path);
qa=IK.data(:,strcmp(IK.colheaders,'ankle_angle_r')).*pi/180;
t_IK=IK.data(:,1);
qa_dot_temp=diff(qa)./diff(t_IK);    qa_dot_temp=[qa_dot_temp; 0];
qa_dot = spline(t_IK,qa_dot_temp,Time);
AnklePower=qa_dot.*Ta;                  % joint power
JointWork=trapz(Time,AnklePower);       % net work
JointPowerPos=AnklePower;     JointPowerPos(AnklePower<0)=0;    PosWork=trapz(Time,JointPowerPos);
JointPowerNeg=AnklePower;     JointPowerNeg(AnklePower>0)=0;    NegWork=trapz(Time,JointPowerNeg);


% To Do: Add metabolic energy computation and add to figure
figure();
% ATendon=[2.5 2.7  2.9 3 4 7 8 9 10 15 35 60 100 150 200];
ATendon=[3 10 20 50];
Cols=hsv(length(ATendon)+10);
for i=1:length(ATendon);
    
    % adjust tendon stiffness
    setup.auxdata.Atendon=ATendon(i);    
    
    % run optimization
    output = gpops2(setup);
    
    % get simulation output
    res=output.result.solution.phase(1);
    Time=res.time;
    MActivation=res.state(:,1:auxdata.NMuscles);
    TForcetilde=res.state(:,auxdata.NMuscles+1:auxdata.NMuscles*2);
    TForce=TForcetilde.*(ones(size(Time))*DatStore.Fiso);
    MExcitation=res.control(:,1:auxdata.NMuscles);
    RActivation=res.control(:,auxdata.NMuscles+1:auxdata.NMuscles+auxdata.Ndof);
    TForce_dottilde=res.control(:,auxdata.NMuscles+auxdata.Ndof+1:end);
    MuscleNames=DatStore.MuscleNames;
    OptInfo=output;
    
    lMTinterp=zeros(length(Time),auxdata.NMuscles);
    VMTinterp=zeros(length(Time),auxdata.NMuscles);
    for m=1:auxdata.NMuscles
        [lMTinterp(:,m),VMTinterp(:,m),~] = SplineEval_ppuval(auxdata.LMTSpline(m),Time,1);
    end
    [lM,lMtilde,vM,vMtilde ] = FiberVelocity_Ftilde(TForcetilde,TForce_dottilde,auxdata.params,lMTinterp,VMTinterp,ATendon(i),S.BoolLinearSpring);
    
    % plot simulation results
    subplot(3,3,1)
    plot(Time,MActivation,'Color',Cols(i,:)); hold on;    
    subplot(3,3,2)
    plot(Time,lMtilde,'Color',Cols(i,:)); hold on;    
    subplot(3,3,3)
    plot(Time,TForce,'Color',Cols(i,:)); hold on;    
    subplot(3,3,4)
    Power=TForce.*vM.*-1;
    plot(Time,Power,'Color',Cols(i,:)); hold on;
    
    PowerNeg=Power; PowerNeg(Power>0)=0;
    PowerPos=Power; PowerPos(Power<0)=0;
    
    % select only stance phase
    i0=find(Time>0.7,1,'first');
    iend=find(Time>1.45,1,'first');
    is=i0:iend;
    
    PosWork=trapz(Time(is),PowerPos(is));
    NegWork=trapz(Time(is),PowerNeg(is));
    
    subplot(3,3,5)
    plot(ATendon(i),PosWork,'*','Color',Cols(i,:)); hold on;
    plot(ATendon(i),NegWork,'*','Color',Cols(i,:)); hold on;
    
    % compute metabolic energy
    FMltilde = get_Flm_tilde(lMtilde,auxdata.Faparam);      % Maarten: not sure if this is the right input
    MuscleMass = 0.5;       % guess of muscle mass (kg)
    TwitchRatio = 50;      % guess of percentage slow twitch fibers  
    [energy_total,energy_am,energy_sl,energy_mech] = ComputeMetabolicEnergy_Umberger2003(MExcitation, MActivation,...
        lMtilde,vMtilde,vM,TForce,MuscleMass,TwitchRatio,10,FMltilde);
    E=trapz(Time(is),energy_total(is).*MuscleMass);
    subplot(3,3,6)
    plot(ATendon(i),E,'*','Color',Cols(i,:)); hold on;
    
    % plot metabolic energy
    subplot(3,3,7)
    plot(Time,energy_total,'Color',Cols(i,:)); hold on;
    
    subplot(3,3,8)
    b=bar(i+4,E); hold on;
    b.FaceColor=Cols(i,:); hold on;
    plot(i+4,PosWork,'*k');
    
    
end

subplot(3,3,1)
xlabel('Time [s]');
ylabel('Muscle activity []');
subplot(3,3,2)
xlabel('Time [s]');
ylabel('LM-tilde []');
subplot(3,3,3)
xlabel('Time [s]');
ylabel('Muscle Force [N]');
set(gca,'YLim',[0 2000]);
subplot(3,3,4)
xlabel('Time [s]');
ylabel('Muscle power [N]');
set(gca,'YLim',[-200 150]);
subplot(3,3,5)
xlabel('ATendon');
ylabel('Muscle fiber work: stance phase');
subplot(3,3,6)
xlabel('ATendon');
ylabel('Metabolic work stance phase ');
subplot(3,3,7)
plot(Time,AnklePower,'k','LineWidth',2);
set(gca,'YLim',[-200 200]);
xlabel('Time [s]');
ylabel('Power [W]');

subplot(3,3,8)
JointWork=trapz(Time(is),AnklePower(is));       % net work
JointPowerPos=AnklePower;     JointPowerPos(AnklePower<0)=0;    PosWork=trapz(Time(is),JointPowerPos(is));
JointPowerNeg=AnklePower;     JointPowerNeg(AnklePower>0)=0;    NegWork=trapz(Time(is),JointPowerNeg(is));
bar(1,JointWork,'k');
bar(2,PosWork,'k');
bar(2,NegWork,'k');
bar(3,abs(PosWork)+abs(NegWork),'k');
ylabel('Joint Work');