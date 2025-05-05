function [tauC,tauK,theta,threshold,R,omega] = taxationLKC_sigma_par(sigma,tauCmax,tauKmax,s)

pi_ew = s.pi_ew;
frac_e = s.frac_e;
skewness = s.skewness;
kurtosis = s.kurtosis;
Ne = s.Ne;

% recalculate model parameters

pi_we = frac_e*pi_ew/(1-frac_e);
Pi_o = [1-pi_we pi_we; pi_ew 1-pi_ew];

if s.Ne == 1 | sigma == 0
    s.Ne = 1;
    s.N = 2;
    PiA = 1;
    logA = 0;
    piA = 1;
else
    % Tanaka-Toda discretization
    cMoments = [0, sigma^2, skewness*sigma^3, kurtosis*sigma^4]; % vector of centered moments
    [logA,piA] = discreteNP(s.Ne,cMoments);
    PiA = repmat(piA,s.Ne,1);
end

s.A = [0 exp(logA)]';
s.Pi = [1-pi_we pi_we*piA; pi_ew*ones(s.Ne,1) (1-pi_ew)*PiA];
s.varpi = [1-frac_e; frac_e*piA'];

% calculate optimal tax rates

eq_sigma = getEq(s);
func = @(tau)(-opttauKobj_L0(tau,tauCmax,eq_sigma.T,s));
options = optimset('TolX',1e-10);
tauK = fminbnd(func,0,tauKmax,options);
[~,s_opt,eq_opt] = opttauKobj_L0(tauK,tauCmax,eq_sigma.T,s);
tauC = s_opt.tau_C;
if sigma == 0
    theta = [eq_opt.theta(1);eq_opt.theta(end)*ones(Ne,1)];
else
    theta = eq_opt.theta;
end
threshold = eq_opt.threshold;
R = eq_opt.R;
omega = eq_opt.omega;

end