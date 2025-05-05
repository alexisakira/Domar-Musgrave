%% optimal labor and capital tax rates with no consumption tax

% maximum tax rates considered

tauLmax = 0.4;
tauKmax = 0.8;

% compute optimal tax rates

func = @(tau)(-opttauKobj(tau,tauLmax,eq.T,s));
options = optimset('TolX',1e-10);
tauK_opt_noC = fminbnd(func,0,tauKmax,options);
[~,s_opt_noC,eq_opt_noC] = opttauKobj(tauK_opt_noC,tauLmax,eq.T,s);

% compute revenue-preserving tax frontier and associated Pareto exponents,
% welfare levels, and aggregate capital and consumption levels

Ntau = 40;
tau_Kvec = linspace(0,tauKmax,Ntau);
tau_Lvec = 0*tau_Kvec;
alphavec = 0*tau_Kvec;
welvec = 0*tau_Kvec;
AggKvec = 0*tau_Kvec;
AggCvec = 0*tau_Kvec;
thetamat = repmat(0*tau_Kvec,s.N,1);
Rvec = 0*tau_Kvec;
omegavec = 0*tau_Kvec;
threshold = 0*tau_Kvec;

s.eps = .0001;

parfor n = 1:Ntau
    [tau_Lvec(n),alphavec(n),welvec(n),AggKvec(n),AggCvec(n),...
        thetamat(:,n),threshold(n),Rvec(n),omegavec(n)] = taxationLK_par(tau_Kvec(n),...
        tauLmax,eq.T,s);
end

% compute truncated wealth distribution at baseline and optimal tax rates

Wgrid_body = (0:.2:10);
Wcdf_body = Wcf2cdf(s,eq,Wgrid_body);
Wcdf_opt_noC_body = Wcf2cdf(s_opt_noC,eq_opt_noC,Wgrid_body);

% compute tail wealth distribution at baseline and optimal tax rates

Wgrid_tail = exp(log(10^1):.1:log(10^3));
Wcdf_tail = Wcf2cdf(s,eq,Wgrid_tail);
Wcdf_opt_noC_tail = Wcf2cdf(s_opt_noC,eq_opt_noC,Wgrid_tail);

% uncomment next line to create backup file
% save('LK_save')

% plot labor income tax rates against capital income tax rates (Figure 4a)

figure(6)
plot(100*tau_Kvec,100*tau_Lvec,'Color',c1)
x0 = xlim;
y0 = ylim;
line(100*[s_opt_noC.tau_K,s_opt_noC.tau_K],[y0(1),100*s_opt_noC.tau_L],...
    'Color',c3,'LineStyle',':')
line([x0(1),100*s_opt_noC.tau_K],100*[s_opt_noC.tau_L,s_opt_noC.tau_L],...
    'Color',c3,'LineStyle',':')
xlabel('Capital income tax rate (\%)')
ylabel('Labor income tax rate (\%)')

% save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_tauK_tauL','-dpdf')

% plot welfare against capital income tax rates (Figure 4b)

figure(7)
plot(100*tau_Kvec,welvec)
x0 = xlim;
y0 = ylim;
line(100*[s_opt_noC.tau_K,s_opt_noC.tau_K],[y0(1),eq_opt_noC.wel_new],...
    'Color',c3,'LineStyle',':')
line([x0(1),100*s_opt_noC.tau_K],[eq_opt_noC.wel_new,...
    eq_opt_noC.wel_new],'Color',c3,'LineStyle',':')
xlabel('Capital income tax rate (\%)')
ylabel('Welfare')

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_tauK_W','-dpdf')

% plot aggregate capital against capital income tax rates (Figure 4c)

figure(8)
plot(100*tau_Kvec,AggKvec,'Color',c1)
x0 = xlim;
y0 = ylim;
AggK_opt = sum(eq_opt_noC.S)-eq_opt_noC.h;
line(100*[s_opt_noC.tau_K,s_opt_noC.tau_K],[y0(1),AggK_opt],...
    'Color',c3,'LineStyle',':')
line([x0(1),100*s_opt_noC.tau_K],[AggK_opt,AggK_opt],'Color',c3,...
    'LineStyle',':')
xlabel('Capital income tax rate (\%)')
ylabel('Capital')

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_tauK_AggK','-dpdf')

% plot aggregate consumption against capital income tax rates (Figure 4d)

figure(9)
plot(100*tau_Kvec,AggCvec,'Color',c1)
x0 = xlim;
y0 = ylim;
AggC_opt = sum(eq_opt_noC.C);
line(100*[s_opt_noC.tau_K,s_opt_noC.tau_K],[y0(1),AggC_opt],...
    'Color',c3,'LineStyle',':')
line([x0(1),100*s_opt_noC.tau_K],[AggC_opt,AggC_opt],'Color',c3,...
    'LineStyle',':')
xlabel('Capital income tax rate (\%)')
ylabel('Consumption')

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_tauK_AggC','-dpdf')

% plot pre-tax interest rate against capital income tax rates (Figure 4e)

figure(103)
yyaxis left
plot(100*tau_Kvec,100*(Rvec-1)./(1-tau_Kvec),'Color',c1)
ylabel('Pre-tax interest rate (\%)','Color','k')
x0 = xlim;
y0 = ylim;
R_opt = 100*(eq_opt_noC.R-1)/(1-s_opt_noC.tau_K);
line(100*[s_opt_noC.tau_K,s_opt_noC.tau_K],[y0(1),R_opt],...
    'Color',c3,'LineStyle',':')
line([x0(1),100*s_opt_noC.tau_K],[R_opt,R_opt],'Color',c3,...
    'LineStyle',':')
yyaxis right
plot(100*tau_Kvec,100*(Rvec-1),'Color',c2,'LineStyle','--')
ylabel('Post-tax interest rate (\%)','Rotation',-90,'Color','k')
x0 = xlim;
y0 = ylim;
R_opt = 100*(eq_opt_noC.R-1);
line(100*[s_opt_noC.tau_K,s_opt_noC.tau_K],[y0(1),R_opt],...
    'Color',c3,'LineStyle',':')
line([100*s_opt_noC.tau_K,x0(2)],[R_opt,R_opt],'Color',c3,...
    'LineStyle',':')
xlabel('Capital income tax rate (\%)')
legend('Pre-tax interest rate','','','Post-tax interest rate','Location','NW')

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_tauK_R','-dpdf')

% plot pre-tax wage against capital income tax rates (Figure 4f)

figure(104)
yyaxis left
plot(100*tau_Kvec,omegavec,'Color',c1)
ylabel('Pre-tax wage','Color','k')
x0 = xlim;
y0 = ylim;
omega_opt = eq_opt_noC.omega;
line(100*[s_opt_noC.tau_K,s_opt_noC.tau_K],[y0(1),omega_opt],...
    'Color',c3,'LineStyle',':')
line([x0(1),100*s_opt_noC.tau_K],[omega_opt,omega_opt],'Color',c3,...
    'LineStyle',':')
yyaxis right
plot(100*tau_Kvec,omegavec.*(1-tau_Lvec),'Color',c2,'LineStyle','--')
ylabel('Post-tax wage','Rotation',-90,'Color','k')
x0 = xlim;
y0 = ylim;
omega_opt = eq_opt_noC.omega*(1-s_opt_noC.tau_L);
line(100*[s_opt_noC.tau_K,s_opt_noC.tau_K],[y0(1),omega_opt],...
    'Color',c3,'LineStyle',':')
line([100*s_opt_noC.tau_K,x0(2)],[omega_opt,omega_opt],'Color',c3,...
    'LineStyle',':')
xlabel('Capital income tax rate (\%)')
legend('Pre-tax wage','','','Post-tax wage','Location','S')

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_tauK_omega','-dpdf')

% plot truncated baseline and optimal wealth distributions (Figure 5a)

figure(10)
plot(Wgrid_body,1-Wcdf_body,'-','Color',c1);
hold on
plot(Wgrid_body,1-Wcdf_opt_noC_body,'--','Color',c2);
xlim([1e-2 1e1])
xlabel('Financial wealth ($w$)')
ylabel('$\mathrm{Pr}(W>w)$')
legend('Baseline','Optimal','Location','NE')

% save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_W_opt_noC_body','-dpdf')

% plot tail baseline and optimal wealth distributions (Figure 5b)

figure(11)
plot(Wgrid_tail,1-Wcdf_tail,'-','Color',c1);
hold on
plot(Wgrid_tail,1-Wcdf_opt_noC_tail,'--','Color',c2);
set(gca,'XScale','log','YScale','log');
xlim([1e1 1e3])
xlabel('Financial wealth ($w$)')
ylabel('$\mathrm{Pr}(W>w)$')
legend('Baseline','Optimal','Location','NE')

% save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_W_opt_noC_tail','-dpdf')

clear tauLmax tauKmax func options tauK_opt_noC Ntau tau_Kvec tau_Lvec...
    alphavec welvec AggKvec x0 y0 fig fig_pos AggW_opt alpha_opt...
    Wgrid_body Wcdf_body Wcdf_opt_noC_body Wgrid_tail Wcdf_tail...
    Wcdf_opt_noC_tail