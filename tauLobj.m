function out = tauLobj(tau,T,s)

s.tau_L = tau;
t = getEq(s);

out = (t.T - T)^2;

end

