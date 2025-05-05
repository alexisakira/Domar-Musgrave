%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MDobj
% (c) 2020 Alexis Akira Toda and Yulong Wang
% 
% Purpose: 
%       Computes objective function of continuously updated minimum distance estimator (CUMDE)
%
% Usage:
%       G = MDobj(xi,p,S)
%
% Inputs:
% xi    - reciprocal of Pareto exponent
% p     - vector of top percentiles
% S     - vector of top income shares
%
% Output:
% G     - objective function
%
% Version 1.0: May 1, 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function G = MDobj(xi,p,S)

%% some error checking
if (xi <= 0)||(xi >= 1)
    error('xi must be between 0 and 1')
    % xi = 1/\alpha must be between 0 and 1
end

if size(p,1) < size(p,2)
    p = p'; % convert p to a column vector
end
if size(S,1) < size(S,2)
    S = S'; % convert S to a column vector
end

K = length(p)-1;

%% compute mu and Sigma

p1 = p(1:end-1); % lower cutoff for percentile group
p2 = p(2:end); % upper cutoff for percentile group

temp1 = (p2.^(1-xi) - p1.^(1-xi))/(1-xi);
mu = temp1;
temp2 = (p2.^(-xi) - p1.^(-xi))/xi;

if xi == 1/2
    temp3 = log(p2./p1);
else
    temp3 = (p2.^(1-2*xi) - p1.^(1-2*xi))/(1-2*xi);
end

sig2 = (2*xi^2/(1-xi))*(temp3 + (p1.^(1-xi)).*temp2 + (2*(p1.^(1-xi)).*(p2.^(1-xi))-p1.^(2-2*xi)-p2.^(2-2*xi))/(2-2*xi));

M = -(xi^2)*temp1*(temp1' + temp2');
U = triu(M,1);
Sigma = diag(sig2) + U + U';

%% compute r, H, Omega

r = mu/mu(end);
r(end) = [];
H = [eye(K-1) -r]/mu(end);
Omega = H*Sigma*(H');

%% compute normalized shares and objective function
S1 = S(1:end-1);
S2 = S(2:end);
temp = S2-S1;
s = temp/temp(end);
s(end) = []; % normalized non-overlapping shares

G = (r-s)'*(Omega\(r-s)); % objective function

end

