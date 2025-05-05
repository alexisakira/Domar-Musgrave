clear all
clc

%% figure formatting

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter', 'latex')
   
set(0,'DefaultTextFontSize', 14)
set(0,'DefaultAxesFontSize', 14)
set(0,'DefaultLineLineWidth',1)

temp = get(gca,'ColorOrder');
c1 = temp(1,:);
c2 = temp(2,:);
c3 = temp(4,:);
c4 = temp(5,:);
clear temp
close all

tic

%% historial tax rates (Figure 1)

historicaltax

%% parameter calibration

s = calibration();

%% compute equilibrium

eq = getEq(s);

%% compute wealth distribution (Figure 2 and Table 2)

wealthdist

%% historical Pareto exponents (Figure 3)

historicalpareto

%% optimal labor and capital taxes with no consumption tax (Figures 4-7)

taxationLK        % Figures 4 and 5

taxationLK_gamma  % Figure 6
taxationLK_sigma  % Figure 7

%% optimal labor, capital and consumption taxes (Figures 8-11)

taxationLKC       % Figures 8 and 9

taxationLKC_gamma % Figure 10
taxationLKC_sigma % Figure 11

%% transition to consumption tax stationary equilibrium (Figure 12)

transition

toc