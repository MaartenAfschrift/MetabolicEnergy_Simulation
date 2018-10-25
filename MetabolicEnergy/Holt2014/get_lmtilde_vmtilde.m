function [lmtilde,vmtilde] = get_lmtilde_vmtilde(strain,time)
%GET_LMTILDE_VMTILDE calculate lmtilde and vmtilde from given strain and lmopt

lmtilde = 1+strain;
pp = spline(time,lmtilde);
dpp = fnder(pp,1);
vmtilde = ppval(dpp,time)./10;

end

