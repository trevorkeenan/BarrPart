function r = gamrnd(a,b,varargin);
%GAMRND Random arrays from gamma distribution.
%   R = GAMRND(A,B) returns an array of random numbers chosen from the
%   gamma distribution with shape parameter A and scale parameter B.  The
%   size of R is the common size of A and B if both are arrays.  If
%   either parameter is a scalar, the size of R is the size of the other
%   parameter.
%
%   R = GAMRND(A,B,M,N,...) or R = GAMRND(A,B,[M,N,...]) returns an
%   M-by-N-by-... array.
%
%   Some references refer to the gamma distribution with a single
%   parameter.  This corresponds to GAMRND with B = 1.
%
%   See also GAMCDF, GAMFIT, GAMINV, GAMLIKE, GAMPDF, GAMSTAT, RANDOM.

%   GAMRND uses a rejection method.

%   References:
%      [1]  Marsaglia, G. and Tsang, W.W. (2000) "A Simple Method for
%           Generating Gamma Variables", ACM Trans. Math. Soft. 26(3):363-372.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2010/10/08 17:23:47 $

if nargin < 2
    error(message('stats:gamrnd:TooFewInputs'));
end

[err, sizeOut] = statsizechk(2,a,b,varargin{:});
if err > 0
    error(message('stats:gamrnd:InputSizeMismatch'));
end

% Return NaN for elements corresponding to illegal parameter values.
a(a < 0) = NaN;
b(b < 0) = NaN;

r = b .* randg(a,sizeOut);
