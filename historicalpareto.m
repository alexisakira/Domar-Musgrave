%% estimate historical wealth Pareto exponents using Toda & Wang (2021)

topProb = [0.01 0.1 0.5 1]/100;
Nc = length(topProb); % number of columns
year = xlsread('.\data\PSZ2020AppendixTablesII(Distrib).xlsx','TE2b','A10:A116');
T = length(year); % number of years
temp = zeros(T,Nc);
temp(:,1) = xlsread('PSZ2020AppendixTablesII(Distrib).xlsx','TE2b','N10:N116'); % top 1%
temp(:,2) = xlsread('PSZ2020AppendixTablesII(Distrib).xlsx','TE2c','B10:B116'); % top 0.5%
temp(:,3) = xlsread('PSZ2020AppendixTablesII(Distrib).xlsx','TE2c','H10:H116'); % top 0.1%
temp(:,4) = xlsread('PSZ2020AppendixTablesII(Distrib).xlsx','TE2c','N10:N116'); % top 0.01%
topShare = fliplr(temp)/100;

alphaHat = zeros(T,1);
for t=1:T
    alphahat = MDestim(topProb,topShare(t,:));
    alphaHat(t) = alphahat;
end

% plot pareto exponents (Figure 3)

figure(4)
plot(year,alphaHat,'-','Color',c1,'LineWidth',1); hold on
plot(year,eq.alpha1*ones(T,1),'--','Color',c2,'LineWidth',1);
text(1920,eq.alpha1,['$\zeta=$' num2str(eq.alpha1,4)],'VerticalAlignment','bottom'); hold off
xlim([year(1) year(end)]);
xlabel('Year')
ylabel('Pareto exponent')
legend('Data','Model')
fontsize(12,"points")

%save figure in pdf format
fig = gcf;
set(fig, 'Units', 'inches', 'Position', [0, 0, 8, 4]); % wider
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_alphaHat','-dpdf')

clear topProb Nc year T temp topShare alphaHat alphahat t fig fig_pos