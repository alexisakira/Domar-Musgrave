% grid of wealth values for Fourier inversion

Wgrid = [linspace(-eq.h+.1,0,1000),exp(log(10^(-2)):.001:log(10^4))];

% grid of wealth values for Pareto extrapolation

Wgrid_tail = exp(log(10^4)+.1:.1:log(10^6));

% compute CDF values

Wcdf = Wcf2cdf(s,eq,Wgrid,Wgrid_tail);

% plot wealth distribution in log-log scale (Figure 2a)

figure(2)
plot(Wgrid,1-Wcdf(1:length(Wgrid)),'-','Color',c1);
hold on
plot(Wgrid_tail,1-Wcdf(length(Wgrid)+1:end),'--','Color',c1);
text(Wgrid(end),1-Wcdf(length(Wgrid)),...
    ['slope $=-\zeta=-$' num2str(eq.alpha1,4)],...
    'HorizontalAlignment','right','VerticalAlignment','top');
set(gca,'XScale','log','YScale','log');
xlim([1e-2 1e6])
xlabel('Financial wealth ($w$)')
ylabel('$\mathrm{Pr}(W>w)$')
legend('Fourier inversion','Pareto extrapolation','Location','NE')

% save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_tailProb','-dpdf')

% compute truncated wealth shares

intQ = trapz(Wgrid,triu(Wcdf(length(Wgrid))-repmat(Wcdf(1:length(...
    Wgrid)),length(Wgrid)-1,1)),2)';
intQ = intQ+Wgrid(1:end-1).*(Wcdf(length(Wgrid))-Wcdf(1:length(Wgrid)-1));
intQ = [intQ,0];

% apply Pareto extrapolation

intQ = intQ+Wgrid(end)*(1-Wcdf(length(Wgrid)))*eq.alpha1/(eq.alpha1-1);

% normalize by aggregate wealth

lorenz = intQ/(sum(eq.S)-eq.h);

% get empirical wealth shares

wsBottom = [-0.2 -0.1 0.2 1.1 2.8 5.6 10.1 17.4]/100; % from Davies et al.
wsTop = [69.2 55.8 33.2 26.5 15.7 7]/100; % from Saez & Zucman
wealthShareHat = flipud([1-wsBottom wsTop]');
topProb = [0.01 0.1 0.5 1 5 10 20 30 40 50 60 70 80 90]/100;

% plot wealth shares in log-log scale (Figure 2b)

figure(3)
plot(1-Wcdf(1:length(Wgrid)),100*lorenz);
hold on
plot(topProb,100*wealthShareHat,'o','Color',c3);
yticks([5 10 20 40 70 100])
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([1e-4 1])
xlabel('Top wealth quantile')
ylabel('Wealth share (\%)')
legend('Model','Data','Location','NW')

% save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_wealthShare','-dpdf')

% get model wealth shares at quantiles in topProb (Table 2)

table2 = [topProb;NaN(1,length(topProb))];
for i = 1:length(topProb)
    n1 = find(Wcdf<=1-topProb(i),1,'last');
    n2 = find(Wcdf>=1-topProb(i),1);
    if n1>=n2
        table2(2,i) = 100*lorenz(round((n1+n2)/2));
    else
        lambda = (log(1-Wcdf(n1))-log(topProb(i)))/(log(1-Wcdf(n1))-...
            log(1-Wcdf(n2)));
        table2(2,i) = 100*exp((1-lambda)*log(lorenz(n1))+...
            lambda*log(lorenz(n2)));
    end
end
save('.\results\table2.mat','table2');

clear fig fig_pos intQ lorenz wsBottom wsTop wealthShareHat topProb i...
    n1 n2 lambda Wgrid Wgrid_tail Wcdf