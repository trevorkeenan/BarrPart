function varargout = dfgetdistributions(varargin)
%DFGETDISTRIBUTIONS A helper function for dffig2m to get structure defining
%the distributions supported by dfittool. 
%   This internal function is used to avoid generating code that that
%   calls dfswitchyard.

%   Copyright 2012 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2012/04/30 03:10:07 $

[varargout{1:nargout}] = dfswitchyard('dfgetdistributions',varargin{:});