function [wel_new,s_opt,eq_opt] = opttauKobj_L0(tau,tauCmax,T,s)

s.tau_K = tau;
s.tau_L = 0;

options = optimset('TolX',1e-10);
func = @(tau)tauCobj(tau,T,s);
tau_C = fminbnd(func,0,tauCmax,options);
s.tau_C = tau_C;

s_opt = s;
eq_opt = getEq(s);
wel_new = eq_opt.wel_new;

end

