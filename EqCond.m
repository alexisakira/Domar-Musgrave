function Out = EqCond(R,omega,s)

% define parameters
beta = s.beta;
gamma = s.gamma;
upsilon = s.upsilon;
tau_K = s.tau_K;
tau_L = s.tau_L;
tau_C = s.tau_C;
alpha = s.alpha;
delta = s.delta;
A = s.A;
Pi = s.Pi;
varpi = s.varpi;
N = s.N;

% solve individual decision problem
r = (1-tau_K)*(A.^(1/alpha)*alpha*(1-alpha)^(1/alpha-1)*omega^(1-1/alpha)-delta);
ell = (A*(1-alpha)/omega).^(1/alpha);
Rk = (1+r)/upsilon; % gross return on capital
Rf = R/upsilon; % effective gross risk-free rate

% compute coefficient of value function and optimal portfolio
[a,theta,~] = get_a(beta,gamma,Pi,Rk,Rf,tau_C,0);

Rtheta = theta*Rk' + Rf*(1-theta)*ones(1,N); % matrix of realized returns
G = beta*Rtheta; % matrix of growth factors
A1 = upsilon*Pi.*G;

if (eigs(A1,1) >= 1)||(R <= upsilon)
    Out = 1e6*[1;1]; % equilibrium condition undefined
    return
end

% now assume A(1) satisfies \rho(A(1)) < 1 and R > \upsilon
h = (1-tau_L)*omega/(1-upsilon/R); % human wealth

excessL = ((1-upsilon)/upsilon)*beta*h*varpi'*((eye(N) - A1)\(theta.*ell)) - 1;
excessB = ((1-upsilon)/upsilon)*beta*varpi'*((eye(N) - A1)\(1-theta)) - 1/R;

Out = [excessL; excessB];

end