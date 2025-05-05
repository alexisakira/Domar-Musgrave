%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MDobj
% (c) 2020 Alexis Akira Toda and Yulong Wang
% 
% Purpose: 
%       Computes continuously updated minimum distance estimator (CUMDE)
%
% Usage:
%       [alpha,G,xi,sigma_alpha,alphaCI] = MDestim(p,S,alpha0,n)
%
% Inputs:
% p     - vector of top percentiles
% S     - vector of top income shares
%
% (Optional inputs):
% alpha0    - initial guess of Pareto exponent
% n         - sample size to compute confidence interval
%
% Outputs:
% alpha         - Pareto exponent
% G             - objective function
% xi            - tail index (reciprocal of Pareto exponent)
% sigma_alpha   - asymptotic standard deviation of alpha
% alphaCI       - 95% confidence interval of alpha
%
% Version 1.0: May 1, 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [alpha,G,xi,sigma_alpha,alphaCI] = MDestim(p,S,alpha0,n)

%% some error checking
if length(p) < 3
    error('length of p must be at least 3')
end
if length(p) ~= length(S)
    error('length of p and S must agree')
end
if any(diff(p) <= 0)
    error('p must be strictly increasing')
end
if (min(p) <= 0)||(max(p) > 1)
    error('elements of p must be between 0 and 1')
end
if any(diff(S) <= 0)
    error('S must be strictly increasing')
end
if (min(S) <= 0)||(max(S) > 1)
    error('elements of S must be between 0 and 1')
end
if size(p,1) < size(p,2)
    p = p'; % convert p to a column vector
end
if size(S,1) < size(S,2)
    S = S'; % convert S to a column vector
end

if (nargin < 3)||isempty(alpha0) % initial guess not provided
    temp = 1-log(S(2)/S(1))/log(p(2)/p(1)); % simple estimator of xi
    if temp > 0 % simple estimator within domain
        xi0 = temp;
    else
        xi0 = 1/2;
    end
else % initial guess provided
    xi0 = 1/alpha0;
end

%% CUMDE

fun = @(xi)MDobj(xi,p,S);
options = optimset('Display','Off');
[xi,G] = fmincon(fun,xi0,[],[],[],[],0.01,0.99,[],options);

alpha = 1/xi; % Pareto exponent
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

%% compute asymptotic variance
temp = (p2.^(1-xi).*log(p2) - p1.^(1-xi).*log(p1))./(p2.^(1-xi) - p1.^(1-xi));
v = temp(end) - temp;
v(end) = [];
R = r.*v;
V = 1/(R'*(Omega\R)); % asymptotic variance
sigma_alpha = sqrt(V)*alpha^2;

%% compute confidence interval if sample size given
if nargout == 5
    if nargin < 4
        error('You need to provide sample size to compute confidence interval')
    end
    alphaCI = [alpha - 1.96*sigma_alpha/sqrt(n), alpha + 1.96*sigma_alpha/sqrt(n)];
end

end

