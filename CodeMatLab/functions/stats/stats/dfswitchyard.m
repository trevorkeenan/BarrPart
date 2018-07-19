function varargout = dfswitchyard(action,varargin)
% DFSWITCHYARD switchyard for Distribution Fitting.
% Helper function for the Distribution Fitting tool

%   Copyright 2003-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.4 $ $Date: 2011/11/09 17:46:37 $

% Calls from Java prefer the if/else version.
% [varargout{1:max(nargout,1)}]=feval(action,varargin{:});
if strcmp(action,'dfhistbins')
    % may appear in generated code from an older release
    action = 'internal.stats.histbins';
end
if nargout==0
	feval(action,varargin{:});
else    
	[varargout{1:nargout}]=feval(action,varargin{:});
end

% The following lines list functions that are called via this
% function from other Statistics Toolbox functions.  These lines
% insure that the compiler will include the functions listed.
%#function mgrp2idx
%#function dfgetdistributions
%#function dfhelpviewer
%#function dfupdatexlim
%#function dfupdateylim
%#function getdsdb
%#function statremovenan
%#function statinsertnan
%#function statkscompute
%#function statparamci
%#function statParallelStore
