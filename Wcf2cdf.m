function cdf = Wcf2cdf(s,eq,Wgrid,Wgrid_tail)
%% convert total log-wealth CF to financial wealth CDF

% get parameters
upsilon = s.upsilon;
varpi = s.varpi;
N = s.N;
Pi = s.Pi;
alpha1 = eq.alpha1;
G = eq.G;
h = eq.h;

% define CF
cf = @(t) (1-upsilon)*(h.^(1i*t)).*reshape(pagemtimes(pagemtimes(...
    varpi',pageinv(eye(N)-upsilon*Pi.*(G.^reshape(1i*t,...
    [1,1,size(t)])))),ones(N,1)),size(t));

% set optional parameters for Fourier inversion
cfoptions.N = 2^18;
cfoptions.T = 10^5;
cfoptions.tolDiff = 1e-6;
cfoptions.isPlot = false;

% apply Fourier inversion over Wgrid
% https://github.com/witkovsky/CharFunTool/blob/master/CF_InvAlgorithms/cf2DistGP.m
[~,cdf,~,~] = cf2DistGP(cf,log(Wgrid+h),[],cfoptions);

% apply Pareto extrapolation over Wgrid_tail
if nargin>3
    cdf_tail = 1-exp(log(1-cdf(end))-alpha1*(log(Wgrid_tail+h)-log(Wgrid(end)+h)));
    cdf = [cdf,cdf_tail];
end

end