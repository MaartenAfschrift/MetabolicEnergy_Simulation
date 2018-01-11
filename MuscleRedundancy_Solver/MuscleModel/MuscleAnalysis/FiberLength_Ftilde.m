% This function computes the muscle fiber length from the normalized tendon
% force

function [lM,lMtilde ] = FiberLength_Ftilde(Ftilde,params,lMT,varargin)

if ~isempty(varargin)
    Atendon=varargin{1};
else
    Atendon=35;
end
lMo = ones(size(Ftilde,1),1)*params(2,:);
lTs = ones(size(Ftilde,1),1)*params(3,:);
alphao = ones(size(Ftilde,1),1)*params(4,:);

% Non-linear tendon
lTtilde = Ftilde./Atendon+1;

% Hill-model relationship
lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);
lMtilde = lM./lMo;
end

