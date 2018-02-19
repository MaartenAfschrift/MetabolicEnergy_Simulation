% This script checks if the smooth approximation gives similar results as
% the non-smooth approximation.

N = 100;
exc = rand(1,N)';
pctst = 50*rand(1,N)';
b = 10000;
f_rec_slow = getPctSTNonSmooth(exc,pctst);
f_rec_slow_sm = getPctSTSmooth(exc,pctst,b);

% CasADi
import casadi.*
exc_SX = SX.sym('exc_SX',N,1);
pctst_SX = SX.sym('pctst_SX',N,1);  
b_SX = SX.sym('b_SX',1); 
rec_slow_sm_SX = getPctSTSmooth(exc_SX,pctst_SX,b_SX);
fgetPctSTSmooth = Function('fgetPctSTSmooth',...
    {exc_SX,pctst_SX,b_SX},{rec_slow_sm_SX});
rec_slow_sm_SX_t = fgetPctSTSmooth(exc,pctst,b);                      
                         
% Check
assert_f_rec_slow = max(abs(f_rec_slow-f_rec_slow_sm));
assert_f_rec_slow_SX = max(abs(f_rec_slow_sm-full(rec_slow_sm_SX_t)));
