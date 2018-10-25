function [t,a] = calculate_activations_excitations(time,u,ta,td,a0)
[t,a] = ode45(@(t,a)activationode(t,a,ta,td,time,u),[0.001 0.4],a0);
