function [y,varargout] = icdf(obj,p,alpha)
%ICDF Inverse cumulative distribution function.
%    Y = ICDF(PD,P) returns an array Y containing the inverse cumulative
%    distribution function (ICDF) for the probability distribution
%    PD, evaluated at values in P.
%
%    See also ProbDist, ProbDistUnivParam, CDF.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $  $Date: 2011/07/20 00:09:08 $

% Check for valid input
if nargin < 2
    error(message('stats:ProbDistUnivParam:icdf:TooFewInputs'));
end

if nargout<=1
   y = paramcall(obj.icdffunc,p,obj.Params);
else
   if nargin<3
       alpha = 0.05;
   end
   [y,varargout{1:nargout-1}] = paramcall(obj.icdffunc,p,obj.Params,obj.ParamCov,alpha);
end