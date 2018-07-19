function varargout = paramcall(fun,x,par,varargin)
%PARAMCALL Helper to call a function with parameters from a vector.

%   Copyright 2008-2010 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2010/10/08 17:28:04 $

% For 1- or 2-parameter calls, do this as efficiently as possible
if isscalar(par)
    [varargout{1:nargout}] = fun(x,par,varargin{:});
elseif numel(par)==2
    [varargout{1:nargout}] = fun(x,par(1),par(2),varargin{:});
else
    pc = num2cell(par);
    [varargout{1:nargout}] = fun(x,pc{:},varargin{:});
end