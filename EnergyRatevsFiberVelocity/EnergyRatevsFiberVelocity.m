% Add necessary folders to path
[parentdir,~,~]=fileparts(pwd);
addpath(genpath(parentdir));

% Load parameters of Muscle-Tendon model
load ActiveFVParameters;
load Faparam;
load PassiveFLParameters;


% -----  SOLEUS  -----%

FMo = 3549;
lMo = 0.05;
lTs = 0.25;
vMmax = 10;

rho = 1059.7;     % Muscle density [kg/m³]
sigma = 0.25e6;   % Specific tension [Pa]
PCSA = FMo/sigma; % Muscle Cross-Sectional Area [m²]
muscle_mass = PCSA*rho*lMo; % Mass [kg]

%% Force Velocity
vM = [0:-0.01:-10]';  % Consider only shortening
vMtilde = vM./vMmax;

e1 = 1.475*ActiveFVParameters(1);
e2 = 0.25*ActiveFVParameters(2);
e3 = ActiveFVParameters(3) + 0.75;
e4 = ActiveFVParameters(4) - 0.027;

FMvtilde = e1*log((e2*vMtilde+e3)+sqrt((e2*vMtilde+e3).^2+1))+e4;
figure(1)
plot(-vMtilde, FMvtilde,'LineWidth',2);
title('Normalised Force-Velocity Multiplier vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Fv/Fvmax');


%% Power Velocity
Power = -1*FMo*FMvtilde.*vM*lMo/muscle_mass;
figure(2)
plot(-vMtilde, Power,'LineWidth',2);
title('Maximal Power vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Power [W/kg]');

%% Energy Rates & Efficiency
% We consider following assumptions
%  -- We have a steady-state muscle activation of 0.5
%  -- The soleus is a slow-twitch muscle consisting for 80% of slow-twitch fibers
%  -- Our muscle is always at optimal fiber length, making the force-length
%     multiplier equal to one and thus making Fce = act*FMo*FMvtilde and Fiso = 1.
%  -- Note that the energy models expect vM to be an absolute speed (need
%     for multiplication with the optimal fiber length compared to formulation in redundancy solver)
exc = 0.5*ones(length(vM),1); act = exc;
pctst = 0.80;
lMtilde = ones(length(vM),1);
Fce = FMo*act.*FMvtilde;
vcemax = vMmax*lMo;
Fiso = ones(length(vM),1);
b = 1000; % smoothness param

%% UMBERGER 2003
[energy_total_Umb2003,energy_am_Umb2003,energy_sl_Umb2003,energy_mech_Umb2003] = ...
    getMetabolicEnergySmooth2003(exc,act,lMtilde,vMtilde,vM*lMo,Fce,muscle_mass, ...
                             pctst,vcemax,Fiso,b);
figure(3)
plot(-vMtilde, energy_total_Umb2003/muscle_mass,-vMtilde,energy_am_Umb2003/muscle_mass,-vMtilde,energy_sl_Umb2003/muscle_mass,-vMtilde,energy_mech_Umb2003/muscle_mass,'LineWidth',2);                         
title('Energy Rates at optimal fiber length and activation 0.5 vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Energy Rate [W/kg]');                  
legend('Total','Activation maintenance heat','Shortening heat','Mechanical')

figure(4)
plot(-vMtilde, energy_mech_Umb2003./energy_total_Umb2003,'LineWidth',2);                         
title('Muscle Efficiency at optimal fiber length and activation 0.5 vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Efficiency [-]');    


% %% UMBERGER 2010 - Identical to 2003 for the shortening range.
% [energy_total_Umb2010,energy_am_Umb2010,energy_sl_Umb2010,energy_mech_Umb2010] = ...
%     getMetabolicEnergySmooth2010(exc,act,lMtilde,vMtilde,vM*lMo,Fce,muscle_mass, ...
%                              pctst,vcemax,Fiso,b);
% figure(5)
% plot(-vMtilde, energy_total_Umb2010/muscle_mass,-vMtilde,energy_am_Umb2010/muscle_mass,-vMtilde,energy_sl_Umb2010/muscle_mass,-vMtilde,energy_mech_Umb2010/muscle_mass,'LineWidth',2);                         
% title('Energy Rates at optimal fiber length and activation 0.5 vs Fiber Velocity')
% xlabel('V/Vmax'); ylabel('Energy Rate [W/kg]');                  
% legend('Total','Activation maintenance heat','Shortening heat','Mechanical')
% 
% figure(6)
% plot(-vMtilde, energy_mech_Umb2010./energy_total_Umb2010,'LineWidth',2);                         
% title('Muscle Efficiency at optimal fiber length and activation 0.5 vs Fiber Velocity')
% xlabel('V/Vmax'); ylabel('Efficiency [-]'); 

% -----  RECTUS FEMORIS  -----%
% ...
% We maintain all muscle parameters equal except the slow twitch percentage
pctst = 0.35;

%% Force Velocity
FMvtilde_FT = e1*log((e2*vMtilde+e3)+sqrt((e2*vMtilde+e3).^2+1))+e4;
figure(5)
plot(-vMtilde, FMvtilde_FT,'LineWidth',2);
title('Normalised Force-Velocity Multiplier vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Fv/Fvmax');     


%% Power Velocity
Power_FT = -1*FMo*FMvtilde_FT.*vM*lMo/muscle_mass;
figure(6)
plot(-vMtilde, Power_FT,'LineWidth',2);
title('Maximal Power vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Power [W/kg]');

%% UMBERGER 2003
Fce = FMo*act.*FMvtilde_FT;
[energy_total_Umb2003_FT,energy_am_Umb2003_FT,energy_sl_Umb2003_FT,energy_mech_Umb2003_FT] = ...
    getMetabolicEnergySmooth2003(exc,act,lMtilde,vMtilde,vM*lMo,Fce,muscle_mass, ...
                             pctst,vcemax,Fiso,b);
figure(7)
plot(-vMtilde, energy_total_Umb2003_FT/muscle_mass,-vMtilde,energy_am_Umb2003_FT/muscle_mass,-vMtilde,energy_sl_Umb2003_FT/muscle_mass,-vMtilde,energy_mech_Umb2003_FT/muscle_mass,'LineWidth',2);                         
title('Energy Rates at optimal fiber length and activation 0.5 vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Energy Rate [W/kg]');                  
legend('Total','Activation maintenance heat','Shortening heat','Mechanical')

figure(8)
plot(-vMtilde, energy_mech_Umb2003_FT./energy_total_Umb2003_FT,'LineWidth',2);                         
title('Muscle Efficiency at optimal fiber length and activation 0.5 vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Efficiency [-]');  


%% Mimick Umberger figures

figure(9)
subplot(2,2,1)
plot(-vMtilde, FMvtilde_FT,-vMtilde, FMvtilde,'LineWidth',2);
title('Normalised Force-Velocity Multiplier vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Fv/Fvmax');    
legend('Fast-Twitch','Slow-Twitch');

subplot(2,2,2)
plot(-vMtilde, Power_FT,-vMtilde, Power,'LineWidth',2);
title('Maximal Power vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Power [W/kg]');
legend('Fast-Twitch','Slow-Twitch');

subplot(2,2,3)
plot(-vMtilde, energy_total_Umb2003_FT/muscle_mass,-vMtilde, energy_total_Umb2003/muscle_mass,'LineWidth',2);                         
title('Energy Rates at optimal fiber length and activation 0.5 vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Energy Rate [W/kg]');                  
legend('Fast-Twitch','Slow-Twitch');


subplot(2,2,4)
plot(-vMtilde, energy_mech_Umb2003_FT./energy_total_Umb2003_FT,-vMtilde, energy_mech_Umb2003./energy_total_Umb2003,'LineWidth',2);                         
title('Muscle Efficiency at optimal fiber length and activation 0.5 vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Efficiency [-]'); 
legend('Fast-Twitch','Slow-Twitch');





% -----  RECTUS FEMORIS  -----%
% ...
% Real rectus femoris
% -----  SOLEUS  -----%

FMo = 1169;
lMo = 0.114;
lTs = 0.31;
vMmax = 10;

rho = 1059.7;     % Muscle density [kg/m³]
sigma = 0.25e6;   % Specific tension [Pa]
PCSA = FMo/sigma; % Muscle Cross-Sectional Area [m²]
muscle_mass = PCSA*rho*lMo; % Mass [kg]





%% Force Velocity
FMvtilde_FT = e1*log((e2*vMtilde+e3)+sqrt((e2*vMtilde+e3).^2+1))+e4;
figure(10)
plot(-vMtilde, FMvtilde_FT,'LineWidth',2);
title('Normalised Force-Velocity Multiplier vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Fv/Fvmax');     


%% Power Velocity
Power_FT = -1*FMo*FMvtilde_FT.*vM*lMo/muscle_mass;
figure(11)
plot(-vMtilde, Power_FT,'LineWidth',2);
title('Maximal Power vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Power [W/kg]');


%% Energy Rates & Efficiency
% We consider following assumptions
%  -- We have a steady-state muscle activation of 0.5
%  -- The soleus is a slow-twitch muscle consisting for 80% of slow-twitch fibers
%  -- Our muscle is always at optimal fiber length, making the force-length
%     multiplier equal to one and thus making Fce = act*FMo*FMvtilde and Fiso = 1.
%  -- Note that the energy models expect vM to be an absolute speed (need
%     for multiplication with the optimal fiber length compared to formulation in redundancy solver)
exc = 0.5*ones(length(vM),1); act = exc;
pctst = 0.35;
lMtilde = ones(length(vM),1);
Fce = FMo*act.*FMvtilde_FT;
vcemax = vMmax*lMo;
Fiso = ones(length(vM),1);
b = 1000; % smoothness param

%% UMBERGER 2003
Fce = FMo*act.*FMvtilde_FT;
[energy_total_Umb2003_FT,energy_am_Umb2003_FT,energy_sl_Umb2003_FT,energy_mech_Umb2003_FT] = ...
    getMetabolicEnergySmooth2003(exc,act,lMtilde,vMtilde,vM*lMo,Fce,muscle_mass, ...
                             pctst,vcemax,Fiso,b);
figure(12)
plot(-vMtilde, energy_total_Umb2003_FT/muscle_mass,-vMtilde,energy_am_Umb2003_FT/muscle_mass,-vMtilde,energy_sl_Umb2003_FT/muscle_mass,-vMtilde,energy_mech_Umb2003_FT/muscle_mass,'LineWidth',2);                         
title('Energy Rates at optimal fiber length and activation 0.5 vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Energy Rate [W/kg]');                  
legend('Total','Activation maintenance heat','Shortening heat','Mechanical')

figure(13)
plot(-vMtilde, energy_mech_Umb2003_FT./energy_total_Umb2003_FT,'LineWidth',2);                         
title('Muscle Efficiency at optimal fiber length and activation 0.5 vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Efficiency [-]');  


%% Mimick Umberger figures

figure(14)
subplot(2,2,1)
plot(-vMtilde, FMvtilde_FT,-vMtilde, FMvtilde,'LineWidth',2);
title('Normalised Force-Velocity Multiplier vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Fv/Fvmax');    
legend('Fast-Twitch','Slow-Twitch');

subplot(2,2,2)
plot(-vMtilde, Power_FT,-vMtilde, Power,'LineWidth',2);
title('Maximal Power vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Power [W/kg]');
legend('Fast-Twitch','Slow-Twitch');

subplot(2,2,3)
plot(-vMtilde, energy_total_Umb2003_FT/muscle_mass,-vMtilde, energy_total_Umb2003/muscle_mass,'LineWidth',2);                         
title('Energy Rates at optimal fiber length and activation 0.5 vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Energy Rate [W/kg]');                  
legend('Fast-Twitch','Slow-Twitch');


subplot(2,2,4)
plot(-vMtilde, energy_mech_Umb2003_FT./energy_total_Umb2003_FT,-vMtilde, energy_mech_Umb2003./energy_total_Umb2003,'LineWidth',2);                         
title('Muscle Efficiency at optimal fiber length and activation 0.5 vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Efficiency [-]'); 
legend('Fast-Twitch','Slow-Twitch');


figure(15)
plot(-vMtilde, energy_total_Umb2003_FT/muscle_mass,-vMtilde,energy_am_Umb2003_FT/muscle_mass,-vMtilde,energy_sl_Umb2003_FT/muscle_mass,-vMtilde,energy_mech_Umb2003_FT/muscle_mass,'LineWidth',2);                         
title('Energy Rates at optimal fiber length and activation 0.5 vs Fiber Velocity')
xlabel('V/Vmax'); ylabel('Energy Rate [W/kg]');                  
legend('Total','Activation maintenance heat','Shortening heat','Mechanical')
