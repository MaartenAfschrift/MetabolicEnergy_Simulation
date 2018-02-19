% This script implement the orderly recruitment model described by Uchida
% et al (2016) and similar to the one proposed by Bhargava et al. (2004)

% Inputs
%   u is muscle excitation
%   pctst is percentage (0,1) slow twitch fibers in muscle (fslow in paper)
% Output
%   f_rec_slow is percentage recruited fibers that are slow-twitch fibers

function f_rec_slow = getPctSTNonSmooth(u,pctst)

u_slow = sin(pi/2*u);
u_fast = 1-cos(pi/2*u);

if u == 0
    f_rec_slow = ones(size(u,1),1);
else
    f_rec_slow = (pctst.*u_slow)./(pctst.*u_slow+(1-pctst).*u_fast);
end

end