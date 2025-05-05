%% sensitivity of optimal taxes (no consumption tax) to risk aversion

% maximum tax rates considered

tauLmax = 0.4;
tauKmax = 0.8;

% range of risk aversion coefficients

gamma = [linspace(1,3.25,23),3.26,linspace(3.27,4.48,10),4.49,linspace(4.5,5,5)];
tauL = 0*gamma;
tauK = tauL;
theta = repmat(tauL,s.N,1);
threshold = tauL;
R = tauL;
omega = tauL;

s.eps = .0001;

% compute optimal taxes

parfor n = 1:length(gamma)
    [tauL(n),tauK(n),theta(:,n),threshold(n),R(n),omega(n)] = taxationLK_gamma_par(gamma(n),tauLmax,tauKmax,s);
end

% uncomment next line to create backup file
% save('LKgamma_save')

% plot optimal tax rates against risk aversion (Figure 6a)

fig = figure(12);
set(fig,'defaultAxesColorOrder',[c1; c2]);
yyaxis left
plot(gamma,100*tauL)
ylabel('Labor income tax rate (\%)','Color','k')
yyaxis right
plot(gamma,100*tauK,'LineStyle','--')
xregion(gamma(find(threshold==1,1)),gamma(find(threshold==1,1,'last')))
ylabel('Capital income tax rate (\%)','Rotation',-90,'Color','k')
xlabel('Risk aversion ($\gamma$)')
legend('Labor income','Capital income','Slackness threshold','Location','W')

%save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_LK_gamma_tax','-dpdf')

% plot equilibrium prices against risk aversion (Figure 6b)

fig = figure(13);
set(fig,'defaultAxesColorOrder',[c1; c2]);
yyaxis left
plot(gamma,100*(R-1))
ylabel('Post-tax interest rate (\%)','Color','k')
yyaxis right
plot(gamma,omega,'LineStyle','--')
xregion(gamma(find(threshold==1,1)),gamma(find(threshold==1,1,'last')))
ylabel('Pre-tax wage','Rotation',-90,'Color','k')
xlabel('Risk aversion ($\gamma$)')
legend('Interest rate','Wage','Slackness threshold','Location','W')

%save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_LK_gamma_prices','-dpdf')

clear tauLmax tauKmax gamma tauL tauK n fig fig_pos