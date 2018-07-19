function x = weibinv(p,a,b)
%WEIBINV Obsolete function
%
%   Use WBLINV in place of WEIBINV.

%Old help text follows.
%
%WEIBINV Inverse of the Weibull cumulative distribution function (cdf).
%   X = WEIBINV(P,A,B) returns the inverse of the Weibull cdf with
%   parameters A and B, at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.    

%   References:
%      [1] J. K. Patel, C. H. Kapadia, and D. B. Owen, "Handbook
%      of Statistical Distributions", Marcel-Dekker, 1976.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2010/10/08 17:27:20 $

if nargin < 3,    
    error(message('stats:weibinv:TooFewInputs')); 
end

[errorcode p a b] = distchck(3,p,a,b);

if errorcode > 0
    error(message('stats:weibinv:InputSizeMismatch'));
end

%   Initialize X to zero.
x = zeros(size(p));

%   Return NaN when the arguments are outside their respective limits.
x(a <= 0 | b <= 0 | p < 0 | p > 1) = NaN;

k = find(a > 0 & b > 0 & p > 0 & p < 1);
if any(k)
    x(k) = (-log(1 - p(k)) ./ a(k)) .^ (1 ./ b(k)); 
end

x(p == 1 & a > 0 & b > 0) = Inf;
