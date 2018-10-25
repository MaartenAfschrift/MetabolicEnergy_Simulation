
function [FMltilde, Fce, Fpe] = get_force_length_velocity_holt(lMtilde, Faparam, Fpparam)

% Active muscle force-length characteristics
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

% force velocity characteristics

vMtilde = vM./vMmax;
e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);

FMvtilde = e1*log((e2*vMtilde+e3)+sqrt((e2*vMtilde+e3).^2+1))+e4;  
Fce = a.*FMltilde.*FMvtilde;

% passive muscle force-length characteristics

e0 = 0.6;
kpe = 4;
t5 = exp(kpe * (lMtilde - 0.10e1) / e0);     %lMtilde
Fpe = ((t5 - 0.10e1) - Fpparam(1)) / Fpparam(2);

end