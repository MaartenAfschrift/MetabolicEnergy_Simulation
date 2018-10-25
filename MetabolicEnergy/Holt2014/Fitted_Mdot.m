function [y] = Fitted_Mdot(x)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

a1 =      0.6854;  
b1 =      0.6314;  
c1 =       1.707;  
a2 =      0.2676;  
b2 =       3.466;  
c2 =      -1.597;  
a3 =      0.1316;  
b3 =       6.151;  
c3 =       -4.57;  
a4 =     0.04492;  
b4 =       8.804;
c4 =      -1.217;
       
y=a1*sin(b1.*x+c1) + a2*sin(b2.*x+c2) + a3*sin(b3.*x+c3) + a4*sin(b4.*x+c4);
                    
end

