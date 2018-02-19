% This script compares the models

N = 100;
act = rand(1,N)';
exc = act;
lMtilde = 2*rand(1,N)';
vMtilde = -1 + 2*rand(1,N)';
vM = -1 + 2*rand(1,N)'*10;
FT = rand(1,N)'*1000;
musclemass = rand(1,N)';
pctst = 50*rand(1,N)';
vcemax = rand(1,N)';
Fiso = rand(1,N)';
b = 100;

[energy_total03,energy_am03,energy_sl03,energy_mech03] = ...
    getMetabolicEnergyNonSmooth2003(exc,act,lMtilde,vMtilde,vM,FT, ...
                                musclemass,pctst,vcemax,Fiso);

[energy_total10,energy_am10,energy_sl10,energy_mech10] = ...
    getMetabolicEnergyNonSmooth2010(exc,act,lMtilde,vMtilde,vM,FT, ...
                                musclemass,pctst,vcemax,Fiso);
                            
                        
figure()
subplot(2,2,1)
plot(energy_total03','k','linewidth',3);
hold on;
plot(energy_total10','r--','linewidth',3);
title('energy: total (w/kg)','Fontsize',16);
set(gca,'Fontsize',16);
subplot(2,2,2)
plot(energy_am03','k','linewidth',3);
hold on;
plot(energy_am10','r--','linewidth',3);
title('energy: activation and maintenance (w/kg)','Fontsize',16);
set(gca,'Fontsize',16);
subplot(2,2,3)
plot(energy_sl03','k','linewidth',3);
hold on;
plot(energy_sl10','r--','linewidth',3);
title('energy: shortening and lengthening (w/kg)','Fontsize',16);
set(gca,'Fontsize',16);
subplot(2,2,4)
plot(energy_mech03','k','linewidth',3);
hold on;
plot(energy_mech10','r--','linewidth',3);
title('energy: mechanical (w/kg)','Fontsize',16);
set(gca,'Fontsize',16);
l = legend('2003','2010');
set(l,'Fontsize',16);