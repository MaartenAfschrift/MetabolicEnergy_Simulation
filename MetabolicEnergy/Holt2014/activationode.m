function adot = activationode(t,a,ta,td,time,u)
b = 0.1;
ui = interp1(time,u,t);
f = 0.5*tanh(b*(ui-a));
d1 = 1/(ta*(0.5+1.5*a));
d2 = (0.5+1.5*a)/td;
adot = (d1.*(f+0.5)+d2.*(-f+0.5)).*(ui-a);