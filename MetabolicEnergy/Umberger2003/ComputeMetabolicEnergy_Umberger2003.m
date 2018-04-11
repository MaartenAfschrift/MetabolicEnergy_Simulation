function [energy_total,energy_am,energy_sl,energy_mech] = ComputeMetabolicEnergy_Umberger2003(exc,act,lmtilde,vtilde,vM,force,musclemass,pctst,vcemax_ft,fiso)
%COMPUTEMETABOLICENERGY_UMBERGER 2003 formules uit Umberger 2003 om het
% metabole energieverbruik op spierniveau te berekenen

% afleiden uit OpenSim
pctft = (100-pctst)/100;
s = 1.5;            % scaling for aerobic activities (1.5) or anaerobic activities (1)

%% stimulations and activations

ind_act_1 = find(exc > act);
ind_act_2 = find(exc <= act);

a(ind_act_1,1) = exc(ind_act_1,1);
a(ind_act_2,1) = (exc(ind_act_2,1)+act(ind_act_2,1))/2;

%% heat activation and maintenanance

ind_lv_1 = find(lmtilde <= 1 & vtilde <= 0);
ind_lv_2 = find(lmtilde <= 1 & vtilde >  0);
ind_lv_3 = find(lmtilde >  1 & vtilde <= 0);
ind_lv_4 = find(lmtilde >  1 & vtilde >  0);

hdotam(ind_lv_1,1) = 1.28*pctft+25;
hdotam(ind_lv_2,1) = 1.28*pctft+25;
hdotam(ind_lv_3,1) = 0.4*(1.28*pctft+25)+0.6.*(1.28*pctft+25).*fiso(ind_lv_3,1);
hdotam(ind_lv_4,1) = 0.4*(1.28*pctft+25)+0.6.*(1.28*pctft+25).*fiso(ind_lv_4,1);

aam = a.^(0.6);
energy_am = hdotam.*aam*s; % activaties die 0 zijn betekenen dat er geen energieverbruik is voor activation and maintenance op dat moment (zie toevoeging paper op 104 L boven)

%% heat shortening and lengthening

coef_hs_ft = 1*153/vcemax_ft;
coef_hs_st = 4*25/(vcemax_ft/2.5); % onduidelijk af te leiden uit de uitleg bij formule (9) en (10) hoe dit af te leiden is uit de constanten AREL en BREL dus gewoon 2.5 genomen die daar ook geciteerd wordt 
coef_hl = 4*coef_hs_st;

hdotsl(ind_lv_1,1) = -1.*coef_hs_st.*vtilde(ind_lv_1,1).*pctst - coef_hs_ft.*vtilde(ind_lv_1,1).*pctft;
hdotsl(ind_lv_2,1) = coef_hl.*vtilde(ind_lv_2,1);
hdotsl(ind_lv_3,1) = -1.*coef_hs_st.*vtilde(ind_lv_3,1).*pctst - coef_hs_ft.*vtilde(ind_lv_3,1).*pctft;
hdotsl(ind_lv_4,1) = coef_hl.*vtilde(ind_lv_4,1);

energy_sl(ind_lv_1,1) = hdotsl(ind_lv_1,1).*(a(ind_lv_1,1).^2).*s;
energy_sl(ind_lv_2,1) = hdotsl(ind_lv_2,1).*a(ind_lv_2,1).*s;
energy_sl(ind_lv_3,1) = hdotsl(ind_lv_3,1).*fiso(ind_lv_3,1).*(a(ind_lv_3,1).^2).*s; % activaties die 0 zijn betekenen dat er geen energieverbruik is voor shortening & lenghtening op dat moment
energy_sl(ind_lv_4,1) = hdotsl(ind_lv_4,1).*fiso(ind_lv_4,1).*a(ind_lv_4,1).*s;

%% mechanical work

energy_mech = -1*force.*vM/musclemass;

%% combine

energy_total = energy_am+energy_sl+energy_mech;     % remark: this is in W/kg => multiple with muscle mass to get actual power


end

