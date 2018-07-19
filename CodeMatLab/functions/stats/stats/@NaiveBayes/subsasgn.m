function [varargout] = subsasgn(varargin)
%SUBSASGN Subscripted reference for a NaiveBayes object.
%   Subscript assignment is not allowed for a NaiveBayes object.

%   Copyright 2008 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2010/10/08 17:27:48 $

error(message('stats:NaiveBayes:subsasgn:NotAllowed'))