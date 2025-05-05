%% transition to new stationary equilibrium with consumption tax

eps = .0001;
Npath = 100;
knots = [1 5 10 20 50 Npath];
knots_add = [12 8 3 6 15 30 2 4 17 40 25 60 70 35 45 80 90 22 13 18 23 11 7 14 9 16 24 19 21];

tau_L = s_opt.tau_L*ones(1,Npath);
tau_K = s_opt.tau_K*ones(1,Npath);
tau_C = s_opt.tau_C*ones(1,Npath);

options = optimset('MaxFunEvals',200000);

for n = 0:length(knots_add)
    if n == 0
        R = linspace(eq.R,eq_opt.R,Npath);
        omega = linspace(eq.omega,eq_opt.omega,Npath);
    else
        knots = sort([knots,knots_add(n)]);
    end
    [x,fval] = fminsearch(@(x)norm(transition_fn(knots,...
        [x(1:length(knots)-1),eq_opt.R],[x(length(knots):end),eq_opt.omega],...
        s,eq,eq_opt,tau_L,tau_K,tau_C,eps)),[R(knots(1:end-1)),...
        omega(knots(1:end-1))],options);
    [exD,theta,S,K,B,h,R,omega,a,r,ell,threshold] = transition_fn(knots,...
        [x(1:length(knots)-1),eq_opt.R],[x(length(knots):end),eq_opt.omega],...
        s,eq,eq_opt,tau_L,tau_K,tau_C,eps);
end

% generate graphs of dynamic responses

time = 1:Npath;
maxT = 50;
time = [0,time(1:maxT)];
R = [eq.R,R(1:maxT)];
omega = [eq.omega,omega(1:maxT)];
theta = [eq.theta,theta(:,1:maxT)];
S = [eq.S,S(:,1:maxT)]; % start-of-period total wealth
K = [s.upsilon*s.Pi'*eq.K,K(:,1:maxT)]; % start-of-period capital
B = [s.upsilon*s.Pi'*eq.B,B(:,1:maxT)]; % start-of-period bond holdings
h = [eq.h,h(1:maxT)];
r = [eq.r,r(:,1:maxT)];
ell = [eq.ell,ell(:,1:maxT)];
tau_L = [s.tau_L,tau_L(1:maxT)];
tau_K = [s.tau_K,tau_K(1:maxT)];
tau_C = [s.tau_C,tau_C(1:maxT)];

% plot of equilibrium interest rate path (Figure 12a)

fig = figure(24);
set(fig,'defaultAxesColorOrder',[c1; c2]);
plot(time,100*(R-1)) % post-tax
ylabel('Post-tax interest rate (\%)')
xlabel('Time (years)')
% save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_trans_R','-dpdf')

% plot of equilibrium wage path (Figure 12b)

fig = figure(25);
set(fig,'defaultAxesColorOrder',[c1; c2]);
plot(time,omega)
ylabel('Pre-tax wage')
xlabel('Time (years)')
% save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_trans_omega','-dpdf')

% end-of-period total wealth

S_end = (1+r).*K+[eq.R,R(1:end-1)].*B+h.*s.varpi;

% consumption

AggC_worker = (1-s.beta)*([1,zeros(1,s.N-1)]*S_end)./(1+tau_C);
AggC_entrepreneur = (1-s.beta)*([0,ones(1,s.N-1)]*S_end)./(1+tau_C);
AggC = AggC_worker+AggC_entrepreneur;

% plot of aggregate consumption by workers and entrepreneurs (Figure 12d)

fig = figure(26);
set(fig,'defaultAxesColorOrder',[c1; c2]);
yyaxis left
plot(time,AggC_worker)
ylabel('Consumption (workers)','Color','k')
yyaxis right
plot(time,AggC_entrepreneur,'LineStyle','--')
ylabel('Consumption (entrepreneurs)','Rotation',-90,'Color','k')
xlabel('Time (years)')
legend('Workers','Entrepreneurs','Location','SE')
% save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_trans_C','-dpdf')

% plot of bond holdings of workers and entrepreneurs (Figure 12e)

fig = figure(27);
set(fig,'defaultAxesColorOrder',[c1; c2; c3]);
plot(time,B(1,:),'-','Color',c1)
hold on
plot(time,sum(B(2:end,:)),'--','Color',c2)
xlabel('Time (years)')
ylabel('Bond holdings')
legend('Workers','Entrepreneurs','Location','best')
% save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_trans_B','-dpdf')

% plot of aggregate capital (Figure 12c)

fig = figure(28);
set(fig,'defaultAxesColorOrder',[c1; c2; c3]);
plot(time,sum(K),'-','Color',c1)
xlabel('Time (years)')
ylabel('Capital')
% save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_trans_K','-dpdf')

% tax revenue by source

revenue_L = tau_L.*omega;
revenue_K = (tau_K./(1-tau_K)).*sum(r.*K);
revenue_C = AggC.*tau_C;
revenue = revenue_L+revenue_K+revenue_C;
revenue_worker = (1-s.frac_e)*revenue_L+AggC_worker.*tau_C;
revenue_entrepreneur = s.frac_e*revenue_L+revenue_K+AggC_entrepreneur.*tau_C;

% plot of tax revenue by source (Figure 12f)

fig = figure(29);
set(fig,'defaultAxesColorOrder',[c1; c2; c3; c4]);
plot(time,revenue_worker,'-','Color',c1)
hold on
plot(time,revenue_entrepreneur,'--','Color',c2)
plot(time,revenue,':','Color',c3)
xlabel('Time (years)')
ylabel('Tax revenue')
legend('Workers','Entrepreneurs','Total','Location','best')
% save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_trans_revenue','-dpdf')

%% political support

Sgrid = [linspace(.1,eq_opt_noC.h,1000),...
    eq_opt_noC.h+exp(log(10^(-2)):.001:log(10^2))]; % grid for total wealth
V_old = zeros(s.N,s.N,length(Sgrid));
V_new = V_old;
cutoff = zeros(s.N,s.N);
vote = cutoff; % matrix to record votes of surviving agents
for J0 = 1:s.N % loop over type in period before regime change
    K0 = (s.beta/s.upsilon)*eq.theta(J0)*Sgrid; % start-of-period capital
    B0 = (s.beta/s.upsilon)*(1-eq.theta(J0))*Sgrid...
        -eq.h/eq.R; % start-of-period bond holdings
    basisvec = eye(s.N);
    basisvec = basisvec(:,J0);
    % Mellin transform of stationary distribution of total wealth
    cf = @(t) (1-s.upsilon)*(eq.h.^(1i*t)).*reshape(pagemtimes(...
        pagemtimes(s.varpi',pageinv(eye(s.N)-s.upsilon*s.Pi.*(eq.G.^...
        reshape(1i*t,[1,1,size(t)])))),basisvec),size(t))/s.varpi(J0);
    for J1 = 1:s.N % loop over type in period of regime change
        S1_old = (1+eq.r(J1))*K0+eq.R*B0+eq.h; % wealth under old regime
        S1_new = (1+r(J1,2))*K0+eq.R*B0+h(2); % wealth under new regime
        V_old(J0,J1,:) = eq.a(J1)*S1_old; % value function under old regime
        V_new(J0,J1,:) = a(J1,2)*S1_new; % value function under new regime
        % uncomment next two lines to plot value function differences
        %figure (s.N*(J0-1)+J1)
        %plot(Sgrid,squeeze(V_new(J0,J1,:)-V_old(J0,J1,:)))
        cutoff(J0,J1) = Sgrid(find(squeeze(V_new(J0,J1,:)-...
            V_old(J0,J1,:))>0,1,'last')); % intersection of value functions
        % Fourier inversion to compute fraction of agents benefitting
        % from regime change
        vote(J0,J1) = Wcf2jointPr(s,eq,J0,cutoff(J0,J1)-eq.h)/s.varpi(J0);
    end
end
vote_total = s.upsilon*sum((s.Pi.*s.varpi).*vote,"all")+...
    (1-s.upsilon)*(a(:,2)*h(2)>eq.a*eq.h)'*s.varpi;
vote_worker = s.upsilon*sum((s.Pi(:,1).*s.varpi).*vote(:,1),"all")...
    /s.varpi(1)+(1-s.upsilon)*(a(1,2)*h(2)>eq.a(1)*eq.h);
vote_entrepreneur = sum((s.Pi(:,2:4).*s.varpi).*vote(:,2:4),"all")...
    /(1-s.varpi(1))+(1-s.upsilon)*(a(2:4,2)*h(2)>eq.a(2:4)*eq.h)'*...
    s.varpi(2:4)/(1-s.varpi(1));

votes = [vote_total vote_worker vote_entrepreneur];
save('.\results\votes.mat','votes');


clear starting_vals knots R_start omega_start Npath time tau_L tau_K...
    tau_C options x fval theta S h R omega a r fig fig_pos AggK AggC...
    AggC_worker AggC_entrepreneur AggB AggB_worker AggB_entrepreneur...
    revenue Sgrid V_old V_new cutoff J0 J1 B0 basisvec cf S1_old...
    S1_new maxT delay K B ell S_end revenue_L revenue_K revenue_C