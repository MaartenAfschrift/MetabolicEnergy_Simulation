function [ fiso,fce ] = get_fce_fiso( a,lMtilde,vMtilde,fmo)
%GET_FCE_FISO the parameters fce and fiso are necessary inputs to compute
% muscle metabolic energy consumption.

%   Fce is defined as the muscle force from the contractile element (length + velocity
%   components but not passive component). Fiso is defined as normalized
%   muscle force from active force-length relationship.


load ActiveFVParameters.mat
Fvparam(1) = 1.475*ActiveFVParameters(1);
Fvparam(2) = 0.25*ActiveFVParameters(2);
Fvparam(3) = ActiveFVParameters(3) + 0.75;
Fvparam(4) = ActiveFVParameters(4) - 0.027;

% Parameters of active muscle force-length characteristic
load Faparam.mat                       

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

fiso = FMtilde1+FMtilde2+FMtilde3;

e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);

FMvtilde = e1*log((e2*vMtilde+e3)+sqrt((e2*vMtilde+e3).^2+1))+e4;

fce = a.*fiso.*FMvtilde.*fmo;

end

