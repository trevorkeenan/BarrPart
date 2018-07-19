function [varargout] = subsasgn(varargin)
%SUBSASGN Subscripted reference for a CLASSREGTREE object.
%   Subscript assignment is not allowed for a CLASSREGTREE object.

%   Copyright 2006-2007 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2010/10/08 17:28:26 $

error(message('stats:classregtree:subsasgn:NotAllowed'))