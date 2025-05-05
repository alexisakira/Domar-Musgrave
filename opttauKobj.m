function [wel_new,s_opt,eq_opt] = opttauKobj(tau,tauLmax,T,s)

s.tau_K = tau;

options = optimset('TolX',1e-10);
func = @(tau)tauLobj(tau,T,s);
tau_L = fminbnd(func,0,tauLmax,options);
s.tau_L = tau_L;

s_opt = s;
eq_opt = getEq(s);
wel_new = eq_opt.wel_new;

end

