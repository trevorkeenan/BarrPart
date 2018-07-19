function classifCheck(C,Sfit,W,cost)

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2011/05/09 01:23:41 $

if ndims(Sfit)>2 || ndims(C)>2
    error(message('stats:classreg:learning:internal:classifCheck:BadDims'));
end
if any(size(Sfit)~=size(C))
    error(message('stats:classreg:learning:internal:classifCheck:SizeScoreClassCountMismatch'));
end

% Check types of Sfit and C
if ~islogical(C) && ~isnumeric(C)
    error(message('stats:classreg:learning:internal:classifCheck:BadC'));
end
if ~isnumeric(Sfit)
    error(message('stats:classreg:learning:internal:classifCheck:BadSfit'));
end

% Check weights
if ~isfloat(W) || ~isvector(W) || length(W)~=size(C,1) || any(W<0)
    error(message('stats:classreg:learning:internal:classifCheck:BadWeights', size( C, 1 )));
end

% Check cost
K = size(C,2);
if ~isempty(cost) && ...
        (ndims(cost)>2 || ~isnumeric(cost) || any(size(cost)~=[K K]) || any(cost(:)<0))
    error(message('stats:classreg:learning:internal:classifCheck:BadCost', K, K));
end
end
