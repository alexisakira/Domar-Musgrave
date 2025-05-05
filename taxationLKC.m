%% optimal labor income, capital income and consumption tax rates

% surface plots use label rotation subroutine available here:
% https://github.com/phymhan/matlab-axis-label-alignment

% maximum tax rates considered

tauLmax = 0.4;
tauKmax = 0.8;
tauCmax = 0.4;

% compute optimal tax rates (no labor income tax)

func = @(tau)(-opttauKobj_L0(tau,tauCmax,eq.T,s));
options = optimset('TolX',1e-10);
tauK_opt = fminbnd(func,0,tauKmax,options);
[~,s_opt,eq_opt] = opttauKobj_L0(tauK_opt,tauCmax,eq.T,s);

% compute frontier of labor and capital income taxes

NtauK = 40;
tau_Kvec = linspace(0,tauKmax,NtauK);
tau_Lvec = 0*tau_Kvec;
parfor n = 1:NtauK
    [tau_Lvec(n),~,~,~] = taxationLK_par(tau_Kvec(NtauK-n+1),tauLmax,...
        eq.T,s);
end

% extend tau_Lvec values back to zero

tau_Lvec = [linspace(0,tau_Lvec(1,1),1+round(tau_Lvec(1,1)/...
    (tau_Lvec(1,2)-tau_Lvec(1,1)))),tau_Lvec(1,2:NtauK)];
NtauL = length(tau_Lvec);

% compute consumption tax, welfare, aggregate capital, aggregate
% consumption, Pareto exponent, leverage, interest rate and wage over grid
% of labor and capital income tax rates

tau_Cmat = NaN(NtauK,NtauL);
welmat = tau_Cmat;
AggKmat = tau_Cmat;
AggCmat = tau_Cmat;
alphamat = tau_Cmat;
thetamat = tau_Cmat;
Rmat = tau_Cmat;
omegamat = tau_Cmat;

parfor n = 1:NtauK
    [tau_Cmat(n,:),welmat(n,:),AggKmat(n,:),AggCmat(n,:),alphamat(n,:),...
        thetamat(n,:),Rmat(n,:),omegamat(n,:)] = taxationLKC_par(...
        tau_Kvec(n),tau_Lvec,tauCmax,n,eq.T,s);
end

% compute curved edges for leverage threshold

edge = [tau_Kvec',NaN(NtauK,6)];
for n = 1:NtauK
    m = find(thetamat(n,:)==1,1);
    if isempty(m)
        edge(n,2) = tau_Lvec(NtauL-n+1);
        edge(n,3) = welmat(n,NtauL-n+1);
        edge(n,4) = AggKmat(n,NtauL-n+1);
        edge(n,5) = AggCmat(n,NtauL-n+1);
        edge(n,6) = Rmat(n,NtauL-n+1);
        edge(n,7) = omegamat(n,NtauL-n+1);
    else
        edge(n,2) = tau_Lvec(m);
        edge(n,3) = welmat(n,m);
        edge(n,4) = AggKmat(n,m);
        edge(n,5) = AggCmat(n,m);
        edge(n,6) = Rmat(n,m);
        edge(n,7) = omegamat(n,m);
    end
end

% compute truncated wealth distribution at baseline and optimal tax rates

Wgrid_body = (0:.2:10);
Wcdf_body = Wcf2cdf(s,eq,Wgrid_body);
Wcdf_opt_body = Wcf2cdf(s_opt,eq_opt,Wgrid_body*(1+s_opt.tau_C));

% compute tail wealth distribution at baseline and optimal tax rates

Wgrid_tail = exp(log(10^1):.1:log(10^3));
Wcdf_tail = Wcf2cdf(s,eq,Wgrid_tail);
Wcdf_opt_tail = Wcf2cdf(s_opt,eq_opt,Wgrid_tail*(1+s_opt.tau_C));

% uncomment next line to create backup file
% save('LKC_save')

% produce surface plot of revenue-preserving tax mixes (Figure 8a)

f = figure(14);
f.Position(3:4) = [450 450];
mesh(100*tau_Kvec,100*tau_Lvec,100*tau_Cmat','FaceAlpha',0)
hold on
xlabel('Cap.\ inc.\ tax (\%)')
ylabel('Labor inc.\ tax (\%)')
zlabel('Cons.\ tax (\%)')
xlim([0 80])
xticks([0 20 40 60 80])
ylim([0 40])
yticks([0 10 20 30 40])
zlim([0 40])
zticks([0 10 20 30 40])
view(160,25)
patchmatx = zeros(2,NtauK-1);
patchmaty = patchmatx;
patchmatz = patchmatx;
for n = 1:NtauK-1;
    patchmatx(:,n) = [100*tau_Kvec(n);100*tau_Kvec(n+1)];
    patchmaty(:,n) = [100*tau_Lvec(end-n+1);100*tau_Lvec(end-n)];
    patchmatz(:,n) = [100*tau_Cmat(n,end-n+1);100*tau_Cmat(n+1,end-n)];
end
patch(patchmatx,patchmaty,patchmatz,patchmatz,'Edgecolor','interp')
scatter3(100*s_opt_noC.tau_K,100*s_opt_noC.tau_L,0,50,'MarkerFaceColor',...
    [0,.75,.75],'MarkerEdgeColor','k','LineWidth',1)
scatter3(100*s_opt.tau_K,0,100*s_opt.tau_C,70,'pentagram','MarkerFaceColor',...
    [.75,0,0],'MarkerEdgeColor','k','LineWidth',1)
scatter3(100*s.tau_K,100*s.tau_L,0,'^','MarkerFaceColor',...
    [.75,0,.75],'MarkerEdgeColor','k','LineWidth',1)
h = rotate3d;
set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
set(gcf, 'ResizeFcn', @align_axislabel)
align_axislabel([], gca)
temp = get(gca,'ylabel');
temp2 = get(temp,'position');
temp2(1) = 0*temp2(1);
set(temp,'position',temp2);
temp = get(gca,'xlabel');
temp2 = get(temp,'position');
temp2 = 1.08*temp2;
set(temp,'position',temp2);
temp = get(gca,'zlabel');
temp2 = get(temp,'position');
temp2(1) = 1.4*temp2(1);
set(temp,'position',temp2);

% save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_KLC','-dpdf')

% produce surface plot of welfare (Figure 8b)

f= figure(15);
f.Position(3:4) = [450 450];
mesh(100*tau_Kvec,100*tau_Lvec,welmat','FaceAlpha',0)
hold on
xlabel('Cap.\ inc.\ tax (\%)')
ylabel('Labor inc.\ tax (\%)')
zlabel('Welfare')
xlim([0 80])
xticks([0 20 40 60 80])
ylim([0 40])
yticks([0 10 20 30 40])
view(160,25)
patchmatx = zeros(2,NtauK-1);
patchmaty = patchmatx;
patchmatz = patchmatx;
for n = 1:NtauK-1
    patchmatx(:,n) = [100*tau_Kvec(n);100*tau_Kvec(n+1)];
    patchmaty(:,n) = [100*edge(n,2);100*edge(n+1,2)];
    patchmatz(:,n) = [edge(n,3);edge(n+1,3)];
end
patch(patchmatx,patchmaty,patchmatz,patchmatz,'Edgecolor','interp')
for n = 1:NtauK-1
    patchmaty(:,n) = [100*tau_Lvec(end-n+1);100*tau_Lvec(end-n)];
    patchmatz(:,n) = [welmat(n,end-n+1);welmat(n+1,end-n)];
end
patch(patchmatx,patchmaty,patchmatz,patchmatz,'Edgecolor','interp')
scatter3(100*s_opt_noC.tau_K,100*s_opt_noC.tau_L,eq_opt_noC.wel_new,50,...
    'MarkerFaceColor',[0,.75,.75],'MarkerEdgeColor','k','LineWidth',1)
scatter3(100*s_opt.tau_K,0,eq_opt.wel_new,70,'pentagram','MarkerFaceColor',...
    [.75,0,0],'MarkerEdgeColor','k','LineWidth',1)
scatter3(100*s.tau_K,100*s.tau_L,eq.wel_new,'^','MarkerFaceColor',...
    [.75,0,.75],'MarkerEdgeColor','k','LineWidth',1)
h = rotate3d;
set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
set(gcf, 'ResizeFcn', @align_axislabel)
align_axislabel([], gca)
temp = get(gca,'ylabel');
temp2 = get(temp,'position');
temp2(1) = 0*temp2(1);
set(temp,'position',temp2);
temp = get(gca,'xlabel');
temp2 = get(temp,'position');
temp2 = 1.08*temp2;
set(temp,'position',temp2);
temp = get(gca,'zlabel');
temp2 = get(temp,'position');
temp2(1) = 1.1*temp2(1);
set(temp,'position',temp2);

% save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_KLWnew','-dpdf')

% produce surface plot of aggregate capital (Figure 8c)

f = figure(16);
f.Position(3:4) = [450 450];
mesh(100*tau_Kvec,100*tau_Lvec,AggKmat','FaceAlpha',0)
hold on
xlabel('Cap.\ inc.\ tax (\%)')
ylabel('Labor inc.\ tax (\%)')
zlabel('Capital')
xlim([0 80])
xticks([0 20 40 60 80])
ylim([0 40])
yticks([0 10 20 30 40])
view(160,25)
patchmatx = zeros(2,NtauK-1);
patchmaty = patchmatx;
patchmatz = patchmatx;
for n = 1:NtauK-1
    patchmatx(:,n) = [100*tau_Kvec(n);100*tau_Kvec(n+1)];
    patchmaty(:,n) = [100*edge(n,2);100*edge(n+1,2)];
    patchmatz(:,n) = [edge(n,4);edge(n+1,4)];
end
patch(patchmatx,patchmaty,patchmatz,patchmatz,'Edgecolor','interp')
for n = 1:NtauK-1
    patchmaty(:,n) = [100*tau_Lvec(end-n+1);100*tau_Lvec(end-n)];
    patchmatz(:,n) = [AggKmat(n,end-n+1);AggKmat(n+1,end-n)];
end
patch(patchmatx,patchmaty,patchmatz,patchmatz,'Edgecolor','interp')
scatter3(100*s_opt_noC.tau_K,100*s_opt_noC.tau_L,sum(eq_opt_noC.S)-eq_opt_noC.h,...
    50,'MarkerFaceColor',[0,.75,.75],'MarkerEdgeColor','k','LineWidth',1)
scatter3(100*s_opt.tau_K,0,sum(eq_opt.S)-eq_opt.h,70,'pentagram',...
    'MarkerFaceColor',[.75,0,0],'MarkerEdgeColor','k','LineWidth',1)
scatter3(100*s.tau_K,100*s.tau_L,sum(eq.S)-eq.h,'^','MarkerFaceColor',...
    [.75,0,.75],'MarkerEdgeColor','k','LineWidth',1)
h = rotate3d;
set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
set(gcf, 'ResizeFcn', @align_axislabel)
align_axislabel([], gca)
temp = get(gca,'ylabel');
temp2 = get(temp,'position');
temp2(1) = 0*temp2(1);
set(temp,'position',temp2);
temp = get(gca,'xlabel');
temp2 = get(temp,'position');
temp2 = 1.08*temp2;
set(temp,'position',temp2);
temp = get(gca,'zlabel');
temp2 = get(temp,'position');
temp2(1) = 1.7*temp2(1);
set(temp,'position',temp2);

% save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_KLAggK','-dpdf')

% produce surface plot of aggregate consumption (Figure 8d)

f = figure(17);
f.Position(3:4) = [450 450];
mesh(100*tau_Kvec,100*tau_Lvec,AggCmat','FaceAlpha',0)
hold on
xlabel('Cap.\ inc.\ tax (\%)')
ylabel('Labor inc.\ tax (\%)')
zlabel('Consumption')
xlim([0 80])
xticks([0 20 40 60 80])
ylim([0 40])
yticks([0 10 20 30 40])
view(160,25)
patchmatx = zeros(2,NtauK-1);
patchmaty = patchmatx;
patchmatz = patchmatx;
for n = 1:NtauK-1
    patchmatx(:,n) = [100*tau_Kvec(n);100*tau_Kvec(n+1)];
    patchmaty(:,n) = [100*edge(n,2);100*edge(n+1,2)];
    patchmatz(:,n) = [edge(n,5);edge(n+1,5)];
end
patch(patchmatx,patchmaty,patchmatz,patchmatz,'Edgecolor','interp')
for n = 1:NtauK-1
    patchmaty(:,n) = [100*tau_Lvec(end-n+1);100*tau_Lvec(end-n)];
    patchmatz(:,n) = [AggCmat(n,end-n+1);AggCmat(n+1,end-n)];
end
patch(patchmatx,patchmaty,patchmatz,patchmatz,'Edgecolor','interp')
scatter3(100*s_opt_noC.tau_K,100*s_opt_noC.tau_L,sum(eq_opt_noC.C),50,...
    'MarkerFaceColor',[0,.75,.75],'MarkerEdgeColor','k','LineWidth',1)
scatter3(100*s_opt.tau_K,0,sum(eq_opt.C),70,'pentagram','MarkerFaceColor',...
    [.75,0,0],'MarkerEdgeColor','k','LineWidth',1)
scatter3(100*s.tau_K,100*s.tau_L,sum(eq.C),'^','MarkerFaceColor',...
    [.75,0,.75],'MarkerEdgeColor','k','LineWidth',1)
h = rotate3d;
set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
set(gcf, 'ResizeFcn', @align_axislabel)
align_axislabel([], gca)
temp = get(gca,'ylabel');
temp2 = get(temp,'position');
temp2(1) = 0*temp2(1);
set(temp,'position',temp2);
temp = get(gca,'xlabel');
temp2 = get(temp,'position');
temp2 = 1.08*temp2;
set(temp,'position',temp2);

% save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_KLAggC','-dpdf')


% produce surface plot of post-tax interest rate (Figure 8e)

f= figure(101);
f.Position(3:4) = [450 450];
mesh(100*tau_Kvec,100*tau_Lvec,100*(Rmat'-1),'FaceAlpha',0)
hold on
xlabel('Cap.\ inc.\ tax (\%)')
ylabel('Labor inc.\ tax (\%)')
zlabel('Interest rate (\%)')
xlim([0 80])
xticks([0 20 40 60 80])
ylim([0 40])
yticks([0 10 20 30 40])
zlim([1.6 1.9])
zticks([1.6 1.7 1.8 1.9])
view(200,25)
patchmatx = zeros(2,NtauK-1);
patchmaty = patchmatx;
patchmatz = patchmatx;
for n = 1:NtauK-1
    patchmatx(:,n) = [100*tau_Kvec(n);100*tau_Kvec(n+1)];
    patchmaty(:,n) = [100*edge(n,2);100*edge(n+1,2)];
    patchmatz(:,n) = [100*(edge(n,6)-1);100*(edge(n+1,6)-1)];
end
patch(patchmatx,patchmaty,patchmatz,patchmatz,'Edgecolor','interp')
for n = 1:NtauK-1
    patchmaty(:,n) = [100*tau_Lvec(end-n+1);100*tau_Lvec(end-n)];
    patchmatz(:,n) = [100*(Rmat(n,end-n+1)-1);100*(Rmat(n+1,end-n)-1)];
end
patch(patchmatx,patchmaty,patchmatz,patchmatz,'Edgecolor','interp')
scatter3(100*s_opt_noC.tau_K,100*s_opt_noC.tau_L,100*(eq_opt_noC.R-1),50,...
    'MarkerFaceColor',[0,.75,.75],'MarkerEdgeColor','k','LineWidth',1)
scatter3(100*s_opt.tau_K,0,100*(eq_opt.R-1),70,'pentagram','MarkerFaceColor',...
    [.75,0,0],'MarkerEdgeColor','k','LineWidth',1)
scatter3(100*s.tau_K,100*s.tau_L,100*(eq.R-1),'^','MarkerFaceColor',...
    [.75,0,.75],'MarkerEdgeColor','k','LineWidth',1)
h = rotate3d;
set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
set(gcf, 'ResizeFcn', @align_axislabel)
align_axislabel([], gca)
temp = get(gca,'ylabel');
temp2 = get(temp,'position');
temp2(1) = 1.05*temp2(1);
temp2(2) = 1*temp2(2);
set(temp,'position',temp2);
temp = get(gca,'xlabel');
temp2 = get(temp,'position');
temp2(1) = 1*temp2(1);
temp2(2) = 1.2*temp2(2);
set(temp,'position',temp2);
temp = get(gca,'zlabel');
temp2 = get(temp,'position');
temp2(1) = 1.1*temp2(1);
set(temp,'position',temp2);

% save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_KLR','-dpdf')


% produce surface plot of pre-tax wage (Figure 8f)

f= figure(102);
f.Position(3:4) = [450 450];
mesh(100*tau_Kvec,100*tau_Lvec,omegamat','FaceAlpha',0)
hold on
xlabel('Cap.\ inc.\ tax (\%)')
ylabel('Labor inc.\ tax (\%)')
zlabel('Wage')
xlim([0 80])
xticks([0 20 40 60 80])
ylim([0 40])
yticks([0 10 20 30 40])
view(160,25)
patchmatx = zeros(2,NtauK-1);
patchmaty = patchmatx;
patchmatz = patchmatx;
for n = 1:NtauK-1
    patchmatx(:,n) = [100*tau_Kvec(n);100*tau_Kvec(n+1)];
    patchmaty(:,n) = [100*edge(n,2);100*edge(n+1,2)];
    patchmatz(:,n) = [edge(n,7);edge(n+1,7)];
end
patch(patchmatx,patchmaty,patchmatz,patchmatz,'Edgecolor','interp')
for n = 1:NtauK-1
    patchmaty(:,n) = [100*tau_Lvec(end-n+1);100*tau_Lvec(end-n)];
    patchmatz(:,n) = [omegamat(n,end-n+1);omegamat(n+1,end-n)];
end
patch(patchmatx,patchmaty,patchmatz,patchmatz,'Edgecolor','interp')
scatter3(100*s_opt_noC.tau_K,100*s_opt_noC.tau_L,eq_opt_noC.omega,50,...
    'MarkerFaceColor',[0,.75,.75],'MarkerEdgeColor','k','LineWidth',1)
scatter3(100*s_opt.tau_K,0,eq_opt.omega,70,'pentagram','MarkerFaceColor',...
    [.75,0,0],'MarkerEdgeColor','k','LineWidth',1)
scatter3(100*s.tau_K,100*s.tau_L,eq.omega,'^','MarkerFaceColor',...
    [.75,0,.75],'MarkerEdgeColor','k','LineWidth',1)
h = rotate3d;
set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
set(gcf, 'ResizeFcn', @align_axislabel)
align_axislabel([], gca)
temp = get(gca,'ylabel');
temp2 = get(temp,'position');
temp2(1) = 0*temp2(1);
set(temp,'position',temp2);
temp = get(gca,'xlabel');
temp2 = get(temp,'position');
temp2 = 1.08*temp2;
set(temp,'position',temp2);
temp = get(gca,'zlabel');
temp2 = get(temp,'position');
temp2(1) = 1.1*temp2(1);
set(temp,'position',temp2);

% save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_KLomega','-dpdf')

% plot truncated baseline and optimal wealth distributions (Figure 9a)

figure(18)
plot(Wgrid_body,1-Wcdf_body,'-','Color',c1);
hold on
plot(Wgrid_body,1-Wcdf_opt_body,'--','Color',c2);
xlim([1e-2 1e1])
xlabel('Effective financial wealth ($w$)')
ylabel('$\mathrm{Pr}(W>w(1+\tau_\mathrm{C}))$')
legend('Baseline','Optimal','Location','NE')

% save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_W_opt_body','-dpdf')

% plot tail baseline and optimal wealth distributions (Figure 9b)

figure(19)
plot(Wgrid_tail,1-Wcdf_tail,'-','Color',c1);
hold on
plot(Wgrid_tail,1-Wcdf_opt_tail,'--','Color',c2);
set(gca,'XScale','log','YScale','log');
xlim([1e1 1e3])
xlabel('Effective financial wealth ($w$)')
ylabel('$\mathrm{Pr}(W>w(1+\tau_\mathrm{C}))$')
legend('Baseline','Optimal','Location','NE')

% save figure in pdf format

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'.\results\fig_W_opt_tail','-dpdf')

clear tauLmax tauKmax tauCmax func options tauK_opt NtauK tau_Kvec...
    tau_Lvec n NtauL tau_Cmat welmat AggKmat alphamat Wgrid_body...
    Wcdf_body Wcdf_opt_body Wgrid_tail Wcdf_tail...
    Wcdf_opt_tail f h temp temp2 fig fig_pos