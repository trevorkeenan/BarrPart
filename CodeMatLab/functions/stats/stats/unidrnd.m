function r = unidrnd(n,varargin)
%UNIDRND Random arrays from the discrete uniform distribution.
%   R = UNIDRND(N) returns an array of random numbers chosen uniformly
%   from the set {1, 2, 3, ... ,N}.  The size of R is the size of N.
%
%   R = UNIDRND(N,MM,NN,...) or R = UNIDRND(N,[MM,NN,...]) returns an
%   MM-by-NN-by-... array.
%
%   See also UNIDCDF, UNIDINV, UNIDPDF, UNIDSTAT, UNIFRND, RANDOM.

%   UNIDRND generates continuous random values, and discretizes them.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 1.1.6.2 $  $Date: 2010/10/08 17:27:03 $

if nargin < 1
    error(message('stats:unidrnd:TooFewInputs')); 
end

[err, sizeOut] = statsizechk(1,n,varargin{:});
if err > 0
    error(message('stats:unidrnd:InputSizeMismatch'));
end

% Return NaN for elements corresponding to illegal parameter values.
n(n <= 0 | round(n) ~= n) = NaN;

r = ceil(n .* rand(sizeOut));
