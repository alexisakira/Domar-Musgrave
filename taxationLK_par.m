function [tau_L,alpha,wel,AggK,AggC,theta,threshold,R,omega] = taxationLK_par(tau_K,tauLmax,T,s)

options = optimset('TolX',1e-10);
s.tau_K = tau_K;
func = @(tau)tauLobj(tau,T,s);
tau_L = fminbnd(func,0,tauLmax,options);
s.tau_L = tau_L;
t = getEq(s);
alpha = t.alpha1;
wel = t.wel_new;
AggK = sum(t.S)-t.h;
AggC = sum(t.C);
theta = t.theta;
threshold = t.threshold;
R = t.R;
omega = t.omega;

end