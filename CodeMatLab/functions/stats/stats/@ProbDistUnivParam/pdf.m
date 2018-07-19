function y=pdf(obj,x)
%PDF Probability density function.
%    Y = PDF(OBJ,X) returns an array Y containing the probability density
%    function (PDF) for the probability distribution object OBJ, evaluated
%    at values in X.
%
%    See also ProbDist, ProbDistUnivParam, PDF.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2011/07/20 00:09:10 $

% Check for valid input
if nargin ~= 2
    error(message('stats:ProbDistUnivParam:pdf:TooFewInputs'));
end

y = paramcall(obj.pdffunc,x,obj.Params);

