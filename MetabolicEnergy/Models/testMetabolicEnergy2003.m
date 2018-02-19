% This script checks if the smooth approximation gives similar results as
% the non-smooth approximation.

N = 100;
act = rand(1,N)';
exc = rand(1,N)';
lMtilde = 2*rand(1,N)';
vMtilde = -1 + 2*rand(1,N)';
vM = -1 + 2*rand(1,N)'*10;
FT = rand(1,N)'*1000;
musclemass = rand(1,N)';
pctst = 50*rand(1,N)';
vcemax = rand(1,N)';
Fiso = rand(1,N)';
b = 10000;

% Non-smooth version
[energy_total,energy_am,energy_sl,energy_mech] = ...
    getMetabolicEnergyNonSmooth2003(exc,act,lMtilde,vMtilde,vM,FT, ...
                                musclemass,pctst,vcemax,Fiso);
                            
% Smooth version
[energy_total_sm,energy_am_sm,energy_sl_sm,energy_mech_sm] = ...
    getMetabolicEnergySmooth2003(exc,act,lMtilde,vMtilde,vM,FT,...
                                 musclemass,pctst,vcemax,Fiso,b);
                         
% CasADi
import casadi.*
act_SX = SX.sym('act_SX',N,1);
exc_SX = SX.sym('exc_SX',N,1);
lMtilde_SX = SX.sym('lMtilde_SX',N,1);
vMtilde_SX = SX.sym('vMtilde_SX',N,1);
vM_SX = SX.sym('vM_SX',N,1);
FT_SX = SX.sym('FT_SX',N,1);
Fiso_SX = SX.sym('Fiso_SX',N,1);
musclemass_SX = SX.sym('musclemass_SX',N,1); 
pctst_SX = SX.sym('pctst_SX',N,1);  
vcemax_SX = SX.sym('vcemax_SX',N,1);
b_SX = SX.sym('b_SX',1); 

[energy_total_sm_SX,energy_am_sm_SX,energy_sl_sm_SX,energy_mech_sm_SX] = ...
    getMetabolicEnergySmooth2003(exc_SX,act_SX,lMtilde_SX,vMtilde_SX,vM_SX,FT_SX,...
                             musclemass_SX,pctst_SX,vcemax_SX,Fiso_SX,b_SX);
fgetMetabolicEnergySmooth2003 = Function('fgetMetabolicEnergySmooth2003',...
    {exc_SX,act_SX,lMtilde_SX,vMtilde_SX,vM_SX,FT_SX,musclemass_SX,pctst_SX,...
    vcemax_SX,Fiso_SX,b_SX},...
    {energy_total_sm_SX,energy_am_sm_SX,energy_sl_sm_SX,energy_mech_sm_SX});

[energy_total_sm_SX_t,energy_am_sm_SX_t,energy_sl_sm_SX_t,energy_mech_sm_SX_t] = ...
    fgetMetabolicEnergySmooth2003(exc,act,lMtilde,vMtilde,vM,FT,musclemass, ...
                             pctst,vcemax,Fiso,b);                      
                         
assertResultsa = max(abs(energy_total-energy_total_sm));
assertResultsb = max(abs(energy_total_sm-full(energy_total_sm_SX_t)));

% figure()
% plot(energy_total);
% hold on;
% plot(energy_total_sm);
% plot(full(energy_total_sm_SX_t));