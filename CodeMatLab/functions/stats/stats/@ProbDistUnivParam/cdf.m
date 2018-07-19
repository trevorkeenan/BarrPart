function [y,varargout] = cdf(obj,x,alpha)
%CDF Cumulative distribution function.
%    Y = CDF(PD,X) returns an array Y containing the cumulative
%    distribution function (CDF) for the probability distribution
%    PD, evaluated at values in X.
%
%    See also ProbDist, ProbDistUnivParam, CDF.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $  $Date: 2011/07/20 00:09:06 $

% Check for valid input
if nargin < 2
    error(message('stats:ProbDistUnivParam:cdf:TooFewInputs'));
end

if nargout<=1
   y = paramcall(obj.cdffunc,x,obj.Params);
else
   if nargin<3
       alpha = 0.05;
   end
   [y,varargout{1:nargout-1}] = paramcall(obj.cdffunc,x,obj.Params,obj.ParamCov,alpha);
end
