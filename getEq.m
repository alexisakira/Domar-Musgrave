function out = getEq(s)
%% compute equilibrium

% get good initial guess by doing grid search
Ng = 10; % number of grid points for initial grid search
lb = [1;1]; % lower bound for R, omega
ub = [1.05;2]; % upper bound for R, omega
Rgrid = linspace(lb(1),ub(1),Ng);
omegagrid = linspace(lb(2),ub(2),Ng);

objMat = zeros(Ng);
for i=1:Ng
    for j=1:Ng
        objMat(i,j) = norm(EqCond(Rgrid(i),omegagrid(j),s));
    end
end

min_objMat = min(objMat(:));
[i,j] = find(objMat == min_objMat);
R0 = Rgrid(i);
omega0 = omegagrid(j);

% compute equilibrium prices starting from initial guess
x0 = [R0;omega0];
options = optimset('TolFun',1e-10,'Display','off');
[x,fval,exitflag] = fsolve(@(x)EqCond(x(1),x(2),s),x0,options);

% repeat with finer grid of start values if equilibrium not found
if exitflag<1
    Ng = 100;
    lb = [1;1]; 
    ub = [1.02;1.5];
    Rgrid = linspace(lb(1),ub(1),Ng);
    omegagrid = linspace(lb(2),ub(2),Ng);
    objMat = zeros(Ng);
    for i=1:Ng
        for j=1:Ng
            objMat(i,j) = norm(EqCond(Rgrid(i),omegagrid(j),s));
        end
    end
    min_objMat = min(objMat(:));
    [i,j] = find(objMat == min_objMat);
    R0 = Rgrid(i);
    omega0 = omegagrid(j);
    x0 = [R0;omega0];
    [x,fval,exitflag] = fsolve(@(x)EqCond(x(1),x(2),s),x0,options);
end

% repeat with even finer grid of start values if equilibrium not found
if exitflag<1
    Ng = 1000;
    lb = [1;1]; 
    ub = [1.02;1.5];
    Rgrid = linspace(lb(1),ub(1),Ng);
    omegagrid = linspace(lb(2),ub(2),Ng);
    objMat = zeros(Ng);
    for i=1:Ng
        for j=1:Ng
            objMat(i,j) = norm(EqCond(Rgrid(i),omegagrid(j),s));
        end
    end
    min_objMat = min(objMat(:));
    [i,j] = find(objMat == min_objMat);
    R0 = Rgrid(i);
    omega0 = omegagrid(j);
    x0 = [R0;omega0];
    [x,fval,exitflag] = fsolve(@(x)EqCond(x(1),x(2),s),x0,options);
end

% repeat with even finer and more localized grid of start values if
% equilibrium not found
if exitflag<1
    lb = [Rgrid(max(i-1,1));omegagrid(max(i-1,1))]; 
    ub = [Rgrid(min(i+1,Ng));omegagrid(min(i+1,Ng))];
    Rgrid = linspace(lb(1),ub(1),Ng);
    omegagrid = linspace(lb(2),ub(2),Ng);
    objMat = zeros(Ng);
    for i=1:Ng
        for j=1:Ng
            objMat(i,j) = norm(EqCond(Rgrid(i),omegagrid(j),s));
        end
    end
    min_objMat = min(objMat(:));
    [i,j] = find(objMat == min_objMat);
    R0 = Rgrid(i);
    omega0 = omegagrid(j);
    x0 = [R0;omega0];
    [x,fval,exitflag] = fsolve(@(x)EqCond(x(1),x(2),s),x0,options);
    if exitflag<1
        disp('Warning: equilibrium not found.');
    end
end

R = x(1);
omega = x(2);

% define parameters
beta = s.beta;
gamma = s.gamma;
upsilon = s.upsilon;
tau_K = s.tau_K;
tau_L = s.tau_L;
tau_C = s.tau_C;
alpha = s.alpha;
delta = s.delta;
eps = s.eps;
A = s.A;
Pi = s.Pi;
varpi = s.varpi;

N = length(A); % number of states

% solve individual decision problem
r = (1-tau_K)*(A.^(1/alpha)*alpha*(1-alpha)^(1/alpha-1)*omega^(1-1/alpha)-delta);
ell = (A*(1-alpha)/omega).^(1/alpha);
Rk = (1+r)/upsilon; % gross return on capital
Rf = R/upsilon; % effective gross risk-free rate

% compute coefficient of value function and optimal portfolio
[a,theta,threshold] = get_a(beta,gamma,Pi,Rk,Rf,tau_C,eps);

Rtheta = theta*Rk' + Rf*(1-theta)*ones(1,N); % matrix of realized returns
G = beta*Rtheta; % matrix of growth factors

f = @(z)(log(eigs(upsilon*Pi.*G.^z,1)));
if f(100) < 0
    alpha1 = Inf;
else
    alpha1 = fzero(f,[1 100]); % upper tail Pareto exponent
end
if f(-100) < 0
    alpha2 = Inf;
else
    alpha2 = -fzero(f,[-100,0]); % lower tail Pareto exponent
end

Aalpha1 = upsilon*Pi.*G.^alpha1;
try
    [v,~] = eigs(Aalpha1',1);
    y = v/sum(v); % top tail type distribution
catch
    y = [];
end

% compute aggregate quantities
h = (1-tau_L)*omega/(1-upsilon/R); % human wealth
A1 = upsilon*Pi.*G;
S = (1-upsilon)*h*((eye(N) - A1')\varpi); % aggregate total wealth vector
C = ((1-beta)/(1+tau_C))*S; % aggregate consumption vector
K = beta*(1/upsilon)*theta.*S; % aggregate capital vector
B = beta*(1/upsilon)*(1-theta).*S-(h/R)*s.varpi; % aggregate bonds vector
T = tau_L*omega+(1-upsilon)*h*((eye(N) - A1')\varpi)'*...
    ((beta*tau_K/(1-tau_K))*((Pi*r).*theta)+...
    ((1-beta)*tau_C/(1+tau_C))*ones(N,1)); % tax revenue
% welfare
if gamma == 1
    wel_new = h*exp(varpi'*log(a));
elseif gamma > 0
    wel_new = h*(varpi'*a.^(1-gamma))^(1/(1-gamma));
else
    error('gamma must be positive')
end

% construct output structure
out.R = R; % gross risk-free rate
out.omega = omega; % wage
out.a = a; % coefficients of value function
out.theta = theta; % optimal portfolio
out.Rtheta = Rtheta; % matrix of gross portfolio return
out.r = r; % vector of returns to physical capital
out.G = G; % matrix of gross growth rates
out.h = h; % human wealth
out.ell = ell; % labor demand
out.alpha1 = alpha1; % upper tail Pareto exponent
out.alpha2 = alpha2; % lower tail Pareto exponent
out.y = y; % upper tail type distribution
out.S = S; % aggregate wealth vector
out.K = K; % aggregate capital vector
out.B = B; % aggregate bonds vector
out.T = T; % tax revenue
out.wel_new = wel_new; % welfare
out.EqError = norm(fval); % equilibrium error
out.aggC = sum(C); % aggregate consumption
out.C = C; % aggregate consumption by ability state
out.threshold = threshold; % flag whether constraint is binding

end
