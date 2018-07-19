function [varargout] = subsasgn(varargin)
%SUBSASGN Subscripted reference for a PROBDIST object.
%   Subscript assignment is not allowed for a PROBDIST object.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2011/05/09 01:27:44 $

error(message('stats:ProbDist:subsasgn:NotAllowed'))