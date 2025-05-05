%% sensitivity of optimal taxes to risk aversion

% maximum tax rates considered

tauCmax = 0.4;
tauKmax = 0.9;

% range of risk aversion coefficients
gamma = [linspace(1,1.48,6),1.49,linspace(1.5,1.74,5),1.75,linspace(1.76,5,28)];
tauC = 0*gamma;
tauK = tauC;
theta = repmat(tauC,s.N,1);
threshold = tauC;
R = tauC;
omega = tauC;

s.eps = .0001; % need smaller threshold than in the other graphs

% compute optimal taxes

parfor n = 1:length(gamma)
    [tauC(n),tauK(n),theta(:,n),threshold(n),R(n),omega(n)] = taxationLKC_gamma_par(gamma(n),tauCmax,tauKmax,s);
end

% uncomment next line to create backup file
% save('LKCgamma_save')

% plot optimal tax rates against risk aversion (Figure 10a)

fig = figure(20);
set(fig,'defaultAxesColorOrder',[c1; c2]);
yyaxis left
plot(gamma,100*tauC)
ylabel('Consumption tax rate (\%)','Color','k')
yyaxis right
plot(gamma,100*tauK,'LineStyle','--')
xregion(gamma(find(threshold==1,1)),gamma(find(threshold==1,1,'last')))
ylabel('Capital income tax rate (\%)','Rotation',-90,'Color','k')
xlabel('Risk aversion ($\gamma$)')
legend('Consumption','Capital income','Slackness threshold','Location','S')

%save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_LKC_gamma_tax','-dpdf')

% plot equilibrium prices against risk aversion (Figure 10b)

fig = figure(21);
set(fig,'defaultAxesColorOrder',[c1; c2]);
yyaxis left
plot(gamma,100*(R-1))
ylabel('Post-tax interest rate (\%)','Color','k')
yyaxis right
plot(gamma,omega,'LineStyle','--')
xregion(gamma(find(threshold==1,1)),gamma(find(threshold==1,1,'last')))
ylabel('Wage','Rotation',-90,'Color','k')
xlabel('Risk aversion ($\gamma$)')
legend('Interest rate','Wage','Slackness threshold','Location','S')

%save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_LKC_gamma_prices','-dpdf')

s.eps = .001;

clear tauCmax tauKmax gamma tauC tauK n fig fig_pos