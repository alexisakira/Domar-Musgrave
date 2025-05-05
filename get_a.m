function [a,theta,threshold] = get_a(beta,gamma,Pi,R,Rf,tau_C,eps)
%% compute coefficient of value function and optimal portfolio

% to get initial guess, conjecture a is multiple of 1
N = length(R);
[~,kappa0,~] = get_theta(gamma,Pi,ones(N,1),R,Rf,eps);
a0 = ((1-beta)/(1+tau_C))^(1-beta)*(beta*kappa0).^beta;

% define fixed point equation; use log transformation to force positivity
func = @(x)(x - log(update_a(beta,gamma,Pi,exp(x),R,Rf,tau_C)));
x0 = log(a0);
options = optimset('Display','off');
x = fsolve(func,x0,options);

a = exp(x);
[theta,~,threshold] = get_theta(gamma,Pi,a,R,Rf,eps);

end