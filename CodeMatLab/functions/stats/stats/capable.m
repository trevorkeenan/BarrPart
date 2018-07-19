function [p, Cp, Cpk] = capable(data,specs)
%CAPABLE Capability indices.
%   CAPABLE(DATA,SPECS) returns the probability that a sample from
%   the process will fall outside the specifications, SPECS.
%   The assumptions are that the measured values in DATA are normally
%   distributed with constant mean and variance and the measurements
%   are statistically independent.
%
%   SPECS is a two element vector containing a lower and upper bound. If
%   there is no lower bound, then use -Inf. Similarly, for no upper bound
%   use Inf.
%
%   [P, CP, CPK] = CAPABLE(DATA,SPECS) also returns the capability 
%   indices CP and CPK. 
%
%   CAPABLE treats NaN values in DATA as missing, and ignores them.
%
%   See also CAPABILITY, HISTFIT.

%   Reference: Montgomery, Douglas, Introduction to Statistical
%   Quality Control, John Wiley & Sons 1991 p. 369-374.

%   Copyright 1993-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.4 $  $Date: 2011/10/15 02:06:11 $

if numel(specs) ~= 2,
   error(message('stats:capable:SpecsBadSize'));
end

lb = specs(1);
ub = specs(2);
if lb > ub
  lb = specs(2);
  ub = specs(1);
end

if lb == -Inf && ub == Inf
   error(message('stats:capable:BadSpecs'));
end

if ~isvector(data)
   error(message('stats:capable:VectorRequired'));
end

data = data(~isnan(data));
mu = mean(data);
sigma = std(data);

if lb == -Inf,
   p = normcdf(-ub,-mu,sigma); % P(t > ub)
   if nargout > 1
      Cp = NaN;
      Cpk = NaN;
      disp(getString(message('stats:capable:text_NoLowerBound')));
   end
elseif ub == Inf,
   p = normcdf(lb,mu,sigma); % P(t < lb)
   if nargout > 1
      Cp = NaN;
      Cpk = NaN;
      disp(getString(message('stats:capable:text_NoUpperBound')));
   end
else  
   p = sum(normcdf([lb -ub],[mu -mu],sigma)); % P(t < lb) + P(t > ub)
   Cp = (ub - lb)./(6.*sigma);
   Cpk = min((ub - mu)./(3.*sigma),(mu-lb)./(3.*sigma));
end


