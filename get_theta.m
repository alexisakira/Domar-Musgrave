function [theta,kappa,threshold] = get_theta(gamma,Pi,a,R,Rf,eps)
%% compute optimal portfolio and kappa

% input arguments
% gamma:    relative risk aversion
% Pi:       transition probability matrix
% a:        coefficients of value functions
% R:        vector of gross return on capital
% Rf:       gross risk-free rate
% eps:      tolerance to determine whether constraint is binding

% output arguments
% theta:        capital share
% kappa:        certainty equivalent
% threshold:    0,1,2 flagging that constraint is binding, on threshhold
%               of binding or slack

N = length(a); % number of states
if size(a,1) < size(a,2)
    error('a must be column vector')
end
if size(R,1) < size(R,2)
    error('R must be column vector')
end

% define CRRA utility
if gamma == 1
    nugam = @(x)(log(x));
    nugaminv = @(y)(exp(y));
elseif gamma > 0
    nugam = @(x)(x.^(1-gamma)/(1-gamma));
    nugaminv = @(y)(((1-gamma)*y).^(1/(1-gamma)));
else
    error('gamma must be positive')
end

theta = zeros(N,1);
kappa = zeros(N,1);

for n = 1:N
    pi = Pi(n,:); % n-th row of Pi
    foc = @(theta)(dot(pi,a.^(1-gamma).*(theta*R + (1-theta)*Rf).^(-gamma).*(R-Rf)));
    if foc(0) <= 0
        thetan = 0;
    elseif foc(1) >= 0
        thetan = 1;
    else
        thetan = fzero(foc,[0,1]);
    end
    if foc(1) > eps
        threshold = 0;
    elseif foc(1) < -eps
        threshold = 2;
    else
        threshold = 1;
    end
    theta(n) = thetan;
    kappa(n) = nugaminv(dot(pi,nugam(a.*(thetan*R + (1-thetan)*Rf))));
end

end

