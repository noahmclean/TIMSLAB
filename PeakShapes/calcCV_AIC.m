function [cv, AIC] = calcCV_AIC(y, B, D, wts, pord, lambda)
%CALCCV returns the cross-validation score for a P-spline
%
%   Based on the CV calculation in psNormal.R from JOPS package in CRAN
%   R authors: Paul Eilers and Brian Marx
%   
%   y are the y-values to be fit (x values are used to make B)
%   B is B-spline basis for P-Spllines
%   D is differencing matrix, of order pord
%   wts is a vector of weights (1/sigma^2)
%   pord is the order of the differencing penalty (usually 2)
%   lambda is the smoothing parameter, and
%   cv returns the cross-validation score

m = length(y); 
n = size(B,2);

P = sqrt(lambda)*D; % for bottom of ls design matrix B+
nix = zeros(n-pord,1); % for bottom of y vector y+

Wplus = diag([wts; nix+1]);
Bplus = [B; P];
yplus = [y; nix];

BTWBinvBTW = (Bplus'*Wplus*Bplus)\(Bplus'*Wplus);
beta = BTWBinvBTW*yplus;
Hplus = Bplus*BTWBinvBTW;

h = diag(Hplus(1:m,1:m));

mu = B*beta;
r = (y - mu)./(1-h); % eqn 3.2, Joy of Splines
cv = sqrt(mean(r.^2)); % eqn 3.1, Joy of Splines

% calculate AIC
Sigma = diag(wts); % covariance matrix
l = -0.5*(log(det(Sigma)) + (y-mu)'*(Sigma\(y-mu)) + m*log(2*pi)); %log-likelihood
ed = sum(h);
AIC = 2*ed - 2*l;



