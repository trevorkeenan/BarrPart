function [m,v] = gamstat(a,b)
%GAMSTAT Mean and variance for the gamma distribution.
%   [M,V] = GAMSTAT(A,B) returns the mean and variance of the gamma
%   distribution with shape and scale parameters A and B, respectively.
%
%   Some references refer to the gamma distribution with a single
%   parameter.  This corresponds to the default of B = 1.
%
%   See also GAMCDF, GAMFIT, GAMINV, GAMLIKE, GAMPDF, GAMRND.

%   References:
%      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
%          Functions, Dover, New York, section 26.1.
%      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 1.1.6.2 $  $Date: 2010/10/08 17:23:48 $

if nargin < 1
    error(message('stats:gamstat:TooFewInputs'));
elseif nargin < 2
    b = 1; 
end

% Return NaN for out of range parameters.
a(a <= 0) = NaN;
b(b <= 0) = NaN;

try
    m = a .* b;
    v = m .* b;
catch
    error(message('stats:gamstat:InputSizeMismatch'));
end
