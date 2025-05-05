%% sensitivity of optimal taxes to volatility

% maximum tax rates considered

tauCmax = 0.4;
tauKmax = 0.9;

% range of volatilities

sigma = [linspace(0,.197,16),.198,linspace(.199,.215,5),.216,linspace(.217,.5,50)];
tauC = 0*sigma;
tauK = tauC;
theta = repmat(tauC,s.N,1);
threshold = tauC;
R = tauC;
omega = tauC;

s.eps = .0001;

% compute optimal taxes

parfor n = 1:length(sigma)
    [tauC(n),tauK(n),theta(:,n),threshold(n),R(n),omega(n)] = taxationLKC_sigma_par(sigma(n),tauCmax,tauKmax,s);
end

% uncomment next line to create backup file
% save('LKCsigma_save')

% plot optimal tax rates against volatility (Figure 11a)

fig = figure(21);
set(fig,'defaultAxesColorOrder',[c1; c2]);
yyaxis left
plot(sigma,100*tauC)
ylabel('Consumption tax rate (\%)','Color','k')
yyaxis right
plot(sigma,100*tauK,'LineStyle','--')
xregion(sigma(find(threshold==1,1)),sigma(find(threshold==1,1,'last')))
ylabel('Capital income tax rate (\%)','Rotation',-90,'Color','k')
xlabel('Volatility of productivity ($\sigma$)')
legend('Consumption','Capital income',['Slackness threshold'],'Location','W')

%save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_LKC_sigma_tax','-dpdf')

% plot equilibrium prices against volatility (Figure 11b)

fig = figure(22);
set(fig,'defaultAxesColorOrder',[c1; c2]);
yyaxis left
plot(sigma,100*(R-1))
ylabel('Post-tax interest rate (\%)','Color','k')
yyaxis right
plot(sigma,omega,'LineStyle','--')
xregion(sigma(find(threshold==1,1)),sigma(find(threshold==1,1,'last')))
ylabel('Wage','Rotation',-90,'Color','k')
xlabel('Volatility of productivity ($\sigma$)')
legend('Interest rate','Wage',['Slackness threshold'],'Location','W')

%save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_LKC_sigma_prices','-dpdf')

clear tauCmax tauKmax sigma tauC tauK n fig fig_pos