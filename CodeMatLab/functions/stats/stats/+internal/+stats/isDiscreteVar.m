function tf = isDiscreteVar(v)
% isDiscreteVar Variable is most naturally treated as discrete.
%    TF = isDiscreteVar(X) returns TRUE if X is categorical, a cell array
%    of strings, logical, or a character array.

%   $Revision: 1.1.6.1 $  $Date: 2011/08/02 14:35:38 $
%   Copyright 2011 The MathWorks, Inc.

tf = (isa(v,'categorical') || iscellstr(v) || islogical(v) || ischar(v));

