function [varargout] = treetest(Tree,varargin)
%TREETEST Obsolete function
%   
%   TREETEST will be removed in a future release. Use CLASSREGTREE/TEST instead.
%
%   See also CLASSREGTREE/TEST.

%   Copyright 1993-2011 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2011/05/09 01:27:06 $

narginchk(1,Inf);

if ~isa(Tree,'classregtree')
    if isa(Tree,'struct')
        Tree = classregtree(Tree);
    else
        error(message('stats:treetest:BadTree'));
    end
end

[varargout{1:max(1,nargout)}] = test(Tree,varargin{:});
