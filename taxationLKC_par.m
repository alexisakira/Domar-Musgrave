function [tau_Cvec,welvec,AggKvec,AggCvec,alphavec,thetavec,Rvec,...
    omegavec] = taxationLKC_par(tau_K,tau_Lvec,tauCmax,m,T,s)

NtauL = length(tau_Lvec);
tau_Cvec = NaN(NtauL,1);
welvec = tau_Cvec;
AggKvec = tau_Cvec;
AggCvec = tau_Cvec;
alphavec = tau_Cvec;
thetavec = tau_Cvec;
Rvec = tau_Cvec;
omegavec = tau_Cvec;
s.tau_K = tau_K;
options = optimset('TolX',1e-10);

for n = 1:NtauL-m+1
    s.tau_L = tau_Lvec(n);
    if n == NtauL-m+1
        tau_Cvec(n) = 0;
        s.tau_C = 0;
    elseif n < NtauL-m+1
        func = @(tau)tauCobj(tau,T,s);
        tau_Cvec(n) = fminbnd(func,0,tauCmax,options);
        s.tau_C = tau_Cvec(n);
    end
    t = getEq(s);
    welvec(n) = t.wel_new;
    AggKvec(n) = sum(t.S)-t.h;
    AggCvec(n) = sum(t.C);
    alphavec(n) = t.alpha1;
    thetavec(n) = t.theta(end);
    Rvec(n) = t.R;
    omegavec(n) = t.omega(end);
end

end