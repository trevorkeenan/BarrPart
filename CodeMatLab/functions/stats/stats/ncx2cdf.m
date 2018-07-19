function p = ncx2cdf(x,v,delta)
%NCX2CDF Non-central chi-square cumulative distribution function (cdf).
%   P = NCX2CDF(X,V,DELTA) Returns the non-central chi-square cdf with V 
%   degrees of freedom and non-centrality parameter, DELTA, at the values 
%   in X.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Some texts refer to this distribution as the generalized Rayleigh,
%   Rayleigh-Rice, or Rice distribution.
%
%   See also NCX2INV, NCX2PDF, NCX2RND, NCX2STAT, CHI2CDF, CDF.

%   Reference:
%      [1]  Evans, Merran, Hastings, Nicholas and Peacock, Brian,
%      "Statistical Distributions, Second Edition", Wiley
%      1993 p. 50-52.
%      [2]  R. Kan and X. Zhao, "Numerical Computation of Noncentral
%      Chi-squared Distribution", unpublished manuscript, 2010.

%   Copyright 1993-2010 The MathWorks, Inc.
%   $Revision: 1.1.8.5 $  $Date: 2010/10/08 17:25:37 $

if nargin < 3
    error(message('stats:ncx2cdf:TooFewInputs')); 
end

[errorcode x v delta] = distchck(3,x,v,delta);

if errorcode > 0
    error(message('stats:ncx2cdf:InputSizeMismatch'));
end

% Initialize P to zero.
if isa(x,'single') || isa(v,'single') || isa(delta,'single')
   p = zeros(size(x),'single');
   eps1 = eps(single(1));
   rmin = realmin('single');
else
   p = zeros(size(x));
   eps1 = eps(1);
   rmin = realmin;
end
k0 = isnan(x) | isnan(v) | isnan(delta);
p(k0) = NaN;
p(x==Inf & ~k0) = 1;

%p(x<=0 & ~k0) = 0;  Not needed as p initialized to 0
p(delta < 0) = NaN;  % can't have negative non-centrality parameter.
p(v < 0) = NaN;      % can't have negative d.f.

% 0 d.f. at x=0
k = v==0 & x==0 & delta>=0 & ~k0;
p(k) = exp(-delta(k)/2);

% Central chi2cdf
k = v>=0 & x>0 & delta==0 & isfinite(x) & ~k0;
p(k) = chi2cdf(x(k),v(k));

% Normal case
todo = find(v>=0 & x>0 & delta>0 & isfinite(x) & ~k0);
delta = delta(todo)/2;
v = v(todo)/2;
x = x(todo)/2;

% Compute Chernoff bounds
e0 = log(rmin);
e1 = log(eps1/4); % 1-eps/4 is the smallest value equal to 1
t = 1 - (v+sqrt(v.^2+4*delta.*x))./(2*x);
q = delta.*t./(1-t) - v.*log(1-t) - t.*x;
peq0 = x<delta+v & q<e0;
peq1 = x>delta+v & q<e1;
%p(todo(peq0)) = 0;  Not needed as p initialized to 0
p(todo(peq1)) = 1;
todo(peq0 | peq1) = [];
x(peq0 | peq1) = [];
v(peq0 | peq1) = [];
delta(peq0 | peq1) = [];

% Find index K of the maximal term in the summation series.
% K1 and K2 are lower and upper bounds for K, respectively.
%
% gammaincratio requires that the 2nd argument be above 1.
%
% Indexing of terms in the summation series starts at 0.
K1 = ceil((sqrt((v+x).^2+4*x.*delta) - (v+x))/2);
K = zeros(size(x));
k1above1 = K1>1;
K2 = floor( delta(k1above1).*gammaincratio(x(k1above1),K1(k1above1)) );
K(k1above1) = K2;

% Find Poisson and Poisson*chi2cdf parts for the maximal terms in the
% summation series.
pois = poisspdf(K,delta);
full = pois.*gammainc(x,v+K);

% Sum the series. First go downward from K and then go upward.
% The term for K is added afterwards - it is not included in either sum.
% Every term is a product of Poisson density and gamma cdf.
sumK = zeros(size(x));

% Downward. poisspdf(k-1,delta)/poisspdf(k,delta) = k/delta
poisterm = pois;
fullterm = full;
keep = K>0 & fullterm>0;
k = K;
while any(keep)
    poisterm(keep) = poisterm(keep).*k(keep)./delta(keep);
    k(keep) = k(keep)-1;
    fullterm(keep) = poisterm(keep).*gammainc(x(keep),v(keep)+k(keep));
    sumK(keep) = sumK(keep) + fullterm(keep);
    keep = keep & k>0 & fullterm>eps(sumK);
end

% Upward. poisspdf(k+1,delta)/poisspdf(k,delta) = delta/(k+1)
poisterm = pois;
fullterm = full;
keep = fullterm>0;
k = K;
while any(keep)
    k(keep) = k(keep)+1;
    poisterm(keep) = poisterm(keep).*delta(keep)./k(keep);
    fullterm(keep) = poisterm(keep).*gammainc(x(keep),v(keep)+k(keep));
    sumK(keep) = sumK(keep) + fullterm(keep);
    keep = keep & fullterm>eps(sumK);
end

% Get probabilities
p(todo) = full + sumK;
p(p>1) = 1;
end

% LocalWords:  RND peq
