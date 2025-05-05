function [tauL,tauK,theta,threshold,R,omega] = taxationLK_gamma_par(gamma,tauLmax,tauKmax,s)

s.gamma = gamma;
eq_gamma = getEq(s);
func = @(tau)(-opttauKobj(tau,tauLmax,eq_gamma.T,s));
options = optimset('TolX',1e-10,'display','off');
[tauK] = fminbnd(func,0,tauKmax,options);
[~,s_opt,eq_opt] = opttauKobj(tauK,tauLmax,eq_gamma.T,s);
tauL = s_opt.tau_L;
theta = eq_opt.theta;
threshold = eq_opt.threshold;
R = eq_opt.R;
omega = eq_opt.omega;

end