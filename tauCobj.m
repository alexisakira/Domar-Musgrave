function out = tauCobj(tau,T,s)

s.tau_C = tau;
t = getEq(s);

out = (t.T - T)^2;

end

