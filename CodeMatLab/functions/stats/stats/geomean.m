function m = geomean(x,dim)
%GEOMEAN Geometric mean.
%   M = GEOMEAN(X) returns the geometric mean of the values in X.  When X
%   is an n element vector, M is the n-th root of the product of the n
%   elements in X.  For a matrix input, M is a row vector containing the
%   geometric mean of each column of X.  For N-D arrays, GEOMEAN operates
%   along the first non-singleton dimension.
%
%   GEOMEAN(X,DIM) takes the geometric mean along dimension DIM of X.
%
%   See also MEAN, HARMMEAN, TRIMMEAN.

%   Copyright 1993-2009 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2010/10/08 17:23:51 $

if any(x(:) < 0)
    error(message('stats:geomean:BadData'))
end

if nargin < 2 || isempty(dim)
    % Figure out which dimension sum will work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

n = size(x,dim);
% Prevent divideByZero warnings for empties, but still return a NaN result.
if n == 0, n = NaN; end

% Take the n-th root of the product of elements of X, along dimension DIM.
if nargin < 2
    m = exp(sum(log(x))./n);
else
    m = exp(sum(log(x),dim)./n);
end
