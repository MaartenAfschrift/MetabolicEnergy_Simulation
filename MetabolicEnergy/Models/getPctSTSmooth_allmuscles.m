% This script implement the orderly recruitment model described by Uchida
% et al (2016) and similar to the one proposed by Bhargava et al. (2004)
% Smooth version with tanh approximation

% Inputs
%   u is muscle excitation
%   pctst is percentage (0,1) slow twitch fibers in muscle (fslow in paper)
%   b is parameter determining transition smoothness
% Output
%   f_rec_slow is percentage recruited fibers that are slow-twitch fibers

function f_rec_slow = getPctSTSmooth_allmuscles(u,pctst,b)

u_slow = sin(pi/2*u);
u_fast = 1-cos(pi/2*u);

u_pos = (0.5 + 0.5*tanh(b*(u)));

f_rec_slow = ones(size(u,1),size(u,2)) + (-ones(size(u,1),size(u,2)) + ...
    (pctst.*u_slow)./(pctst.*u_slow+(1-pctst).*u_fast)).*u_pos;

end