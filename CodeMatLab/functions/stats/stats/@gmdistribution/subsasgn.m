function [varargout] = subsasgn(varargin)
%SUBSASGN Subscripted reference for a gmdistribution object.
%   Subscript assignment is not allowed for a gmdistribution object.

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2010/10/08 17:28:43 $

error(message('stats:gmdistribution:subsasgn:NotAllowed'))