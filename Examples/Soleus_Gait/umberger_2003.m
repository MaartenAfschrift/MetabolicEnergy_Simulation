clear all
close all
clc

%% input

DO = load('C:\Users\u0098084\Documents\Projecten\ATSECW\Processing\Christianne_Platteau\20171113\DynamicOptimization\CP_20171113_walking_6MWD_08_DO.mat');
modelpath = 'C:\Users\u0098084\Documents\Projecten\ATSECW\Processing\Christianne_Platteau\20171113\Models\gait2392_simbody_ATSECW_CP_slacklength.osim';
muscle = 'soleus_r';
act = DO.MActivation.stride1(2:end,34);
lmtilde = DO.lMtilde.stride1(2:end,34);
vtilde = diff(DO.lMtilde.stride1(:,34))./diff(DO.Time.stride1);
% figure()
% plot(vtilde); hold on
time = DO.Time.stride1(2:end);
[b,a] = butter(4,10/((length(time)/(time(end)-time(1)))/2), 'low');
vtilde = filtfilt(b,a,vtilde);
% plot(vtilde);
force = DO.TForce.stride1(2:end,34);


fiso = []; % de kracht-lengte karakteristiek uit de statische optimalisatie (uit DO) genomen en geplakt, daarna spline naar de tijdsframes van de SO
% fiso = spline(linspace(1,106,106), fiso, linspace(1,106,315))';
s = 1.5; % scaling for aerobic activities (1.5) or anaerobic activities (1)

% afleiden uit OpenSim
[musclemass, pctst, vcemax_ft] = umberger_2003_inputopensim(modelpath, muscle);
pctft = 1-pctst;

%% heat activation and maintenanance

ind_lv_1 = find(lmtilde <= 1 & vtilde <= 0);
ind_lv_2 = find(lmtilde <= 1 & vtilde >  0);
ind_lv_3 = find(lmtilde >  1 & vtilde <= 0);
ind_lv_4 = find(lmtilde >  1 & vtilde >  0);

hdotam(ind_lv_1,1) = 1.28*pctft+25;
hdotam(ind_lv_2,1) = 1.28*pctft+25;
hdotam(ind_lv_3,1) = 0.4*(1.28*pctft+25)+0.6.*(1.28*pctft+25).*fiso(ind_lv_3,1);
hdotam(ind_lv_4,1) = 0.4*(1.28*pctft+25)+0.6.*(1.28*pctft+25).*fiso(ind_lv_4,1);

aam = act.^(0.6);

energy_am = hdotam.*aam*s; % activaties die 0 zijn betekenen dat er geen energieverbruik is voor activation and maintenance op dat moment (zie toevoeging paper op 104 L boven)

%% heat shortening and lengthening

coef_hs_ft = 1*153/vcemax_ft;
coef_hs_st = 4*25/(vcemax_ft/2.5); % onduidelijk af te leiden uit de uitleg bij formule (9) en (10) hoe dit af te leiden is uit de constanten AREL en BREL dus gewoon 2.5 genomen die daar ook geciteerd wordt 
coef_hl = 4*coef_hs_st;

hdotsl(ind_lv_1,1) = -1.*coef_hs_st.*vtilde(ind_lv_1,1).*pctst - coef_hs_ft.*vtilde(ind_lv_1,1).*pctft;
hdotsl(ind_lv_2,1) = coef_hl.*vtilde(ind_lv_2,1);
hdotsl(ind_lv_3,1) = -1.*coef_hs_st.*vtilde(ind_lv_3,1).*pctst - coef_hs_ft.*vtilde(ind_lv_3,1).*pctft;
hdotsl(ind_lv_4,1) = coef_hl.*vtilde(ind_lv_4,1);

energy_sl(ind_lv_1,1) = hdotsl(ind_lv_1,1).*(act(ind_lv_1,1).^2).*s;
energy_sl(ind_lv_2,1) = hdotsl(ind_lv_2,1).*act(ind_lv_2,1).*s;
energy_sl(ind_lv_3,1) = hdotsl(ind_lv_3,1).*fiso(ind_lv_3,1).*(act(ind_lv_3,1).^2).*s; % activaties die 0 zijn betekenen dat er geen energieverbruik is voor shortening & lenghtening op dat moment
energy_sl(ind_lv_4,1) = hdotsl(ind_lv_4,1).*fiso(ind_lv_4,1).*act(ind_lv_4,1).*s;

%% mechanical work

energy_mech = -1*force.*vtilde/musclemass;

%% combine

energy_total = energy_am+energy_sl+energy_mech;
