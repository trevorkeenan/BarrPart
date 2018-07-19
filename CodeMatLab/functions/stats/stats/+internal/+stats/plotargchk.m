function plotargchk(varargin)
% Internal utility to check plotting arguments

%   Copyright 2012 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2012/05/08 20:43:13 $

if nargin==0
    return
end

if mod(nargin,2)==1
    m = message('stats:internal:parseArgs:WrongNumberArgs');
    throwAsCaller(MException(m.Identifier, '%s', getString(m)));
end

if ~iscellstr(varargin(1:2:end))
    m = message('stats:internal:parseArgs:IllegalParamName');
    throwAsCaller(MException(m.Identifier, '%s', getString(m)));
end
