%% sensitivity of optimal taxes (no consumption tax) to volatility

% maximum tax rates considered

tauLmax = 0.4;
tauKmax = 0.8;

% range of volatilities

sigma = [linspace(0,.252,25),.253,linspace(.254,.274,5),.275,linspace(.276,.4,10)];
tauL = 0*sigma;
tauK = tauL;
theta = repmat(tauL,s.N,1);
threshold = tauL;
R = tauL;
omega = tauL;

s.eps = .0001;

% compute optimal taxes

parfor n = 1:length(sigma)
    [tauL(n),tauK(n),theta(:,n),threshold(n),R(n),omega(n)] = taxationLK_sigma_par(sigma(n),tauLmax,tauKmax,s);
end

% uncomment next line to create backup file
% save('LKsigma_save')

% plot optimal tax rates against volatility (Figure 7a)

fig = figure(13);
set(fig,'defaultAxesColorOrder',[c1; c2]);
yyaxis left
plot(sigma,100*tauL)
ylabel('Labor income tax rate (\%)','Color','k')
yyaxis right
plot(sigma,100*tauK,'LineStyle','--')
xregion(sigma(find(threshold==1,1)),sigma(find(threshold==1,1,'last')))
ylabel('Capital income tax rate (\%)','Rotation',-90,'Color','k')
xlabel('Volatility of productivity ($\sigma$)')
legend('Labor income','Capital income','Slackness threshold','Location','W')

%save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_LK_sigma_tax','-dpdf')

% plot equilibrium prices against volatility (Figure 7b)

fig = figure(14);
set(fig,'defaultAxesColorOrder',[c1; c2]);
yyaxis left
plot(sigma,100*(R-1))
ylabel('Post-tax interest rate (\%)','Color','k')
yyaxis right
plot(sigma,omega,'LineStyle','--')
xregion(sigma(find(threshold==1,1)),sigma(find(threshold==1,1,'last')))
ylabel('Pre-tax wage','Rotation',-90,'Color','k')
xlabel('Volatility of productivity ($\sigma$)')
legend('Interest rate','Wage','Slackness threshold','Location','W')

%save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_LK_sigma_prices','-dpdf')

clear tauLmax tauKmax sigma tauL tauK n fig fig_pos