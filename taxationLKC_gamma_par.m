function [tauC,tauK,theta,threshold,R,omega] = taxationLKC_gamma_par(gamma,tauCmax,tauKmax,s)

s.gamma = gamma;
eq_gamma = getEq(s);
func = @(tau)(-opttauKobj_L0(tau,tauCmax,eq_gamma.T,s));
options = optimset('TolX',1e-10);
tauK = fminbnd(func,0,tauKmax,options);
[~,s_opt,eq_opt] = opttauKobj_L0(tauK,tauCmax,eq_gamma.T,s);
tauC = s_opt.tau_C;
theta = eq_opt.theta;
threshold = eq_opt.threshold;
R = eq_opt.R;
omega = eq_opt.omega;

end