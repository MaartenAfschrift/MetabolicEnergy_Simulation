
function [FMltilde] = get_Flm_tilde_Holt(lMtilde,Faparam)

% Active muscle force-length characterist
b11 = Faparam(1);
b21 = Faparam(2);
b31 = Faparam(3);
b41 = Faparam(4);
b12 = Faparam(5);
b22 = Faparam(6);
b32 = Faparam(7);
b42 = Faparam(8);

b13 = 0.1;
b23 = 1;
b33 = 0.5*sqrt(0.5);
b43 = 0;
num3 = lMtilde-b23;
den3 = b33+b43*lMtilde;
FMtilde3 = b13*exp(-0.5*num3.^2./den3.^2);

num1 = lMtilde-b21;
den1 = b31+b41*lMtilde;
FMtilde1 = b11*exp(-0.5*num1.^2./den1.^2);

num2 = lMtilde-b22;
den2 = b32+b42*lMtilde;
FMtilde2 = b12*exp(-0.5*num2.^2./den2.^2);

FMltilde = FMtilde1+FMtilde2+FMtilde3;


end