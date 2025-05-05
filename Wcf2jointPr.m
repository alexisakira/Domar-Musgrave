function jointPr = Wcf2jointPr(s,eq,n,Wgrid,Wgrid_tail)
%% convert total log-wealth conditional CF to Pr(W<=w,J=n)

% get parameters
upsilon = s.upsilon;
varpi = s.varpi;
N = s.N;
Pi = s.Pi;
alpha1 = eq.alpha1;
G = eq.G;
h = eq.h;

% conditional CF of total log-wealth
basis_n = eye(N);
basis_n = basis_n(:,n);
cf = @(t) (1/varpi(n))*(1-upsilon)*(h.^(1i*t)).*reshape(pagemtimes(pagemtimes(...
    varpi',pageinv(eye(N)-upsilon*Pi.*(G.^reshape(1i*t,...
    [1,1,size(t)])))),basis_n),size(t));

% set optional parameters for Fourier inversion
cfoptions.N = 2^18;
cfoptions.T = 10^5;
cfoptions.tolDiff = 1e-6;
cfoptions.isPlot = false;

% apply Fourier inversion over Wgrid
% https://github.com/witkovsky/CharFunTool/blob/master/CF_InvAlgorithms/cf2DistGP.m
[~,cdf,~,~] = cf2DistGP(cf,log(Wgrid+h),[],cfoptions);

% apply Pareto extrapolation over Wgrid_tail
if nargin>4
    cdf_tail = 1-exp(log(1-cdf(end))-alpha1*(log(Wgrid_tail+h)-log(Wgrid(end)+h)));
    cdf = [cdf,cdf_tail];
end

% multiply by state probability to obtain joint probability

jointPr = varpi(n)*cdf;

end