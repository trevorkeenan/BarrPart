function y=pdf(obj,x)
%PDF Probability density function.
%    Y = PDF(PD,X) returns an array Y containing the probability density
%    function (PDF) for the probability distribution PD, evaluated at
%    values in X.
%
%    See also PROBDISTUNIVKERNEL, PDF.

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2010/10/08 17:28:00 $

% Check for valid input
if nargin ~= 2
    error(message('stats:ProbDistUnivKernel:pdf:TooFewInputs'));
end

ksinfo = obj.ksinfo;
y = dfswitchyard('statkscompute','pdf',x,true,length(x),obj.BandWidth, ...
                 ksinfo.L,ksinfo.U,ksinfo.weight,[],obj.Kernel,...
                 ksinfo.ty,[],ksinfo.foldpoint,ksinfo.maxp);

