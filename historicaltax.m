% historical capital income tax rates

data = xlsread('.\data\taxrates2023.xlsx');

year = data(:,1);
t_cap = data(:,2);
t_corp = data(:,3);
t_K = data(:,4);

% plot capital income tax rates (Figure 1)

figure(1)
stairs(year,t_cap,'--','Color',c1,'LineWidth',1);
hold on
stairs(year,t_corp,':','Color',c3,'LineWidth',1);
stairs(year,t_K,'-','Color',c2,'LineWidth',1);
xlim([min(year) max(year)])
ylim([0 100])
xlabel('Year')
ylabel('Top marginal tax rate (\%)')
legend('Capital gains','Corporate income','Effective capital income')
fontsize(12,"points")

%save figure in pdf format
fig = gcf;
set(fig, 'Units', 'inches', 'Position', [0, 0, 8, 4]); % wider
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_taxrates','-dpdf')

clear data year t_cap t_corp t_K fig fig_pos