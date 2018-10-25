% This script checks if the smooth approximation gives similar results as
% the non-smooth approximation.

N = 100;
act = rand(1,N)';
exc = rand(1,N)';
lMtilde = 2*rand(1,N)';
vM = 4.*rand(1,N)' - 2;
Fce = rand(1,N)'*10;
Fpass = rand(1,N)'*10;
musclemass = rand(1,N)';
pctst = 50*rand(1,N)';
vcemax = rand(1,N)';
Fiso = rand(1,N)';
Fmax = rand(1,N)'*1000;
b = 10000;
modelmass = 100;

% Non-smooth version
[energy_total,Adot,Mdot,Sdot,Wdot,energy_model] = ...
    getMetabolicEnergyNonSmooth2004(exc,act,lMtilde,vM,Fce,Fpass, ...
                                musclemass,pctst,Fiso,Fmax,modelmass);
                            
% Smooth version
[energy_total_sm,Adot_sm,Mdot_sm,Sdot_sm,Wdot_sm,energy_model_sm] = ...
    getMetabolicEnergySmooth2004(exc,act,lMtilde,vM,Fce,Fpass, ...
                                musclemass,pctst,Fiso,Fmax,modelmass,b);
                         
% CasADi version
import casadi.*
act_SX = SX.sym('act_SX',N,1);
exc_SX = SX.sym('exc_SX',N,1);
lMtilde_SX = SX.sym('lMtilde_SX',N,1);
vM_SX = SX.sym('vM_SX',N,1);
Fce_SX = SX.sym('Fce_SX',N,1);
Fpass_SX = SX.sym('Fce_SX',N,1);
Fiso_SX = SX.sym('Fiso_SX',N,1);
musclemass_SX = SX.sym('musclemass_SX',N,1); 
pctst_SX = SX.sym('pctst_SX',N,1);  
Fmax_SX = SX.sym('Fmax_SX',N,1);
b_SX = SX.sym('b_SX',1); 

[energy_total_sm_SX,Adot_sm_SX,Mdot_sm_SX,Sdot_sm_SX,Wdot_sm_SX,...
    energy_mod_sm_SX] = getMetabolicEnergySmooth2004(exc_SX,act_SX,...
    lMtilde_SX,vM_SX,Fce_SX,Fpass_SX,musclemass_SX,pctst_SX,Fiso_SX,...
    Fmax_SX,modelmass,b_SX);
fgetMetabolicEnergySmooth2004=Function('fgetMetabolicEnergySmooth2004',...
    {exc_SX,act_SX,lMtilde_SX,vM_SX,Fce_SX,Fpass_SX,musclemass_SX,...
    pctst_SX,Fiso_SX,Fmax_SX,b_SX},{energy_total_sm_SX,Adot_sm_SX,...
    Mdot_sm_SX,Sdot_sm_SX,Wdot_sm_SX,energy_mod_sm_SX});

[energy_total_sm_SX_t,Adot_sm_SX_t,Mdot_sm_SX_t,Sdot_sm_SX_t,...
    Wdot_sm_SX_t,energy_mod_sm_SX_t] = fgetMetabolicEnergySmooth2004(...
    exc,act,lMtilde,vM,Fce,Fpass,musclemass,pctst,Fiso,Fmax,b);                      
                         
assertResultsa = max(abs(energy_total-energy_total_sm));
assertResultsb = max(abs(energy_total_sm-full(energy_total_sm_SX_t)));
