function s = calibration()
%% define parameter values of the model

beta = .96;      % discount factor
gamma = 3;        % relative risk aversion
upsilon = 1-1/40; % survival probability

tau_K = 1-(1-.238)*(1-.21);      % capital income tax rate
tau_L = .248;     % labor income tax rate
tau_C = 0;        % consumption tax rate

alpha = .36;     % capital share
delta = .08;     % capital depreciation

pi_ew = .0192;     % firm exit rate
frac_e = .115;   % fraction of entrepreneurs

Ne = 5; % number of entrepreneurial states

eps = .001; % tolerance for determining whether constraint is binding

% calibrate parameters from Fagereng et al. (Table 3, risky financial
% asset)
sigma = .2473; % standard deviation
skewness = -.08; % skewness
kurtosis = 6.22; % kurtosis

pi_we = frac_e*pi_ew/(1-frac_e);
Pi_o = [1-pi_we pi_we; pi_ew 1-pi_ew]; % transition probability matrix for occupation

if Ne == 1
    PiA = 1;
    logA = 0;
    piA = 1;
else
    % Tanaka-Toda discretization
    cMoments = [0, sigma^2, skewness*sigma^3, kurtosis*sigma^4]; % vector of centered moments
    [logA,piA] = discreteNP(Ne,cMoments);
    PiA = repmat(piA,Ne,1);
end

% define transition probability matrix and productivities

A = [0 exp(logA)]'; % vector of TFP
Pi = [1-pi_we pi_we*piA; pi_ew*ones(Ne,1) (1-pi_ew)*PiA];
varpi = [1-frac_e; frac_e*piA'];
N = length(A); % number of states

% define structure
s.beta = beta;
s.gamma = gamma;
s.upsilon = upsilon;
s.tau_K = tau_K;
s.tau_L = tau_L;
s.tau_C = tau_C;
s.alpha = alpha;
s.delta = delta;
s.A = A;
s.Pi = Pi;
s.varpi = varpi;
s.N = N;
s.pi_ew = pi_ew;
s.frac_e = frac_e;
s.Ne = Ne;
s.sigma = sigma;
s.skewness = skewness;
s.kurtosis = kurtosis;
s.eps = eps;