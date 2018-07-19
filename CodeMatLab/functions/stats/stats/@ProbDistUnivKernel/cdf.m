function y=cdf(obj,x)
%CDF Cumulative distribution function.
%    Y = CDF(PD,X) returns an array Y containing the cumulative
%    distribution function (CDF) for the probability distribution
%    PD, evaluated at values in X.
%
%    See also PROBDISTUNIVKERNEL, CDF.

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2010/10/08 17:27:58 $

% Check for valid input
if nargin ~= 2
    error(message('stats:ProbDistUnivKernel:cdf:TooFewInputs'));
end

ksinfo = obj.ksinfo;
y = dfswitchyard('statkscompute','cdf',x,true,length(x),obj.BandWidth, ...
                 ksinfo.L,ksinfo.U,ksinfo.weight,[],obj.Kernel,...
                 ksinfo.ty,[],ksinfo.foldpoint,ksinfo.maxp);
