%% convert price paths to excess demand paths

function [exD,theta,S,K,B,h,R,omega,a,r,ell,threshold] = transition_fn(...
    knots,R_knots,omega_knots,s,eq,eq_new,tau_L,tau_K,tau_C,eps)

% interpolate price paths

Npath = knots(end);
R = spline(knots,R_knots,1:Npath);
omega = spline(knots,omega_knots,1:Npath);

% compute returns on capital and optimal labor demand

r = (1-tau_K).*((s.A.^(1/s.alpha))*s.alpha*((1-s.alpha)^...
    (1/s.alpha-1))*(omega.^(1-1/s.alpha))-s.delta);

ell = (s.A*(1-s.alpha)*(omega.^(-1))).^(1/s.alpha);

% backward recursion to compute paths of a, theta and h

a = zeros(s.N,Npath);
a(:,end) = eq_new.a;

theta = zeros(s.N,Npath);
theta(:,end) = eq_new.theta;

h = zeros(1,Npath);
h(end) = eq_new.h;

threshold = zeros(1,Npath);
threshold(end) = 2;

for n = Npath-1:-1:1
    [theta(:,n),kappa,threshold(n)] = get_theta(s.gamma,s.Pi,a(:,n+1),...
        (1+r(:,n+1))/s.upsilon,R(n)/s.upsilon,eps);
    a(:,n) = (((1-s.beta)/(1+tau_C(n)))^(1-s.beta))*...
        ((s.beta*kappa).^s.beta);
    h(n) = (s.upsilon/R(n))*h(n+1)+(1-tau_L(n))*omega(n);
end

% forward recursion to compute excess demand paths

K = zeros(s.N,Npath);
B = zeros(s.N,Npath);
S = zeros(s.N,Npath);
K(:,1) = s.upsilon*s.Pi'*eq.K;
B(:,1) = s.upsilon*s.Pi'*eq.B;
S(:,1) = (1+r(:,1)).*K(:,1)+eq.R*B(:,1)+h(1)*s.varpi;

exB = zeros(1,Npath);
exL = zeros(1,Npath);

for n = 2:Npath
    K(:,n) = s.beta*s.Pi'*(theta(:,n-1).*S(:,n-1));
    B(:,n) = s.Pi'*(s.beta*((1-theta(:,n-1)).*S(:,n-1))-...
        s.upsilon*(h(n)/R(n-1))*s.varpi);
    S(:,n) = (1+r(:,n)).*K(:,n)+R(n-1)*B(:,n)+h(n)*s.varpi;
    exB(n-1) = sum((s.beta/s.upsilon)*((1-theta(:,n-1)).*S(:,n-1))-...
        (h(n)/R(n-1))*s.varpi);
    exL(n-1) = (s.beta/s.upsilon)*(theta(:,n-1).*ell(:,n-1))'*S(:,n-1)-1;
end

exD = [exB,exL];

end