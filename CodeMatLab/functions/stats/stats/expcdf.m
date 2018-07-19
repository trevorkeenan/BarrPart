function [p,plo,pup] = expcdf(x,mu,pcov,alpha)
%EXPCDF Exponential cumulative distribution function.
%   P = EXPCDF(X,MU) returns the cdf of the exponential distribution with
%   mean parameter MU, evaluated at the values in X.  The size of P is
%   the common size of the input arguments.  A scalar input functions as a
%   constant matrix of the same size as the other input.
%
%   The default value for MU is 1.
%
%   [P,PLO,PUP] = EXPCDF(X,MU,PCOV,ALPHA) produces confidence bounds
%   for P when the input parameter MU is an estimate.  PCOV is the
%   variance of the estimated MU.  ALPHA has a default value of 0.05, and
%   specifies 100*(1-ALPHA)% confidence bounds.  PLO and PUP are arrays of
%   the same size as P containing the lower and upper confidence bounds.
%   The bounds are based on a normal approximation for the distribution of
%   the log of the estimate of MU.  You can get a more accurate set of
%   bounds simply by using EXPFIT to get a confidence interval for MU,
%   and evaluating EXPCDF at the lower and upper end points of that interval.
%
%   See also CDF, EXPFIT, EXPINV, EXPLIKE, EXPPDF, EXPRND, EXPSTAT.

%   References:
%     [1] Lawless, J.F. (1982) Statistical Models and Methods for Lifetime Data, Wiley,
%         New York.
%     [2} Meeker, W.Q. and L.A. Escobar (1998) Statistical Methods for Reliability Data,
%         Wiley, New York.
%     [3] Crowder, M.J., A.C. Kimber, R.L. Smith, and T.J. Sweeting (1991) Statistical
%         Analysis of Reliability Data, Chapman and Hall, London.

%     Copyright 1993-2011 The MathWorks, Inc. 
%     $Revision: 1.1.8.4 $  $Date: 2011/05/09 01:24:54 $

if nargin<1
    error(message('stats:expcdf:TooFewInputsX'));
end
if nargin < 2
    mu = 1;
end

% More checking if we need to compute confidence bounds.
if nargout>1
   if nargin<3
      error(message('stats:expcdf:TooFewInputsVariance'));
   end
   if numel(pcov)~=1
      error(message('stats:expcdf:BadVarianceScalar'));
   end
   if nargin<4
      alpha = 0.05;
   elseif ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
      error(message('stats:expcdf:BadAlpha'));
   end
end

% Return NaN for out of range parameters.
mu(mu <= 0) = NaN;

try
    z = x./mu;
catch
    error(message('stats:expcdf:InputSizeMismatch'));
end

% Force a zero for negative x.
z(z < 0) = 0;
p = 1 - exp(-z);

% Compute confidence bounds if requested.
if nargout>=2
   % Work on log scale.
   logz = log(z);
   if pcov<0
      error(message('stats:expcdf:BadVarianceNonNeg'));
   end
   normz = -norminv(alpha/2);
   halfwidth = normz * sqrt(pcov ./ (mu.^2));
   zlo = logz - halfwidth;
   zup = logz + halfwidth;
   
   % Convert back to original scale
   plo = 1 - exp(-exp(zlo));
   pup = 1 - exp(-exp(zup));
end
