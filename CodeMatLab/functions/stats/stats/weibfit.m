function [phat, pci] = weibfit(data,alpha)
%WEIBFIT Obsolete function
%
%   Use WBLFIT in place of WEIBFIT.

%Old help text follows.
%
%WEIBFIT Parameter estimates and confidence intervals for Weibull data.
%   WEIBFIT(DATA,ALPHA) Returns the maximum likelihood estimates of the  
%   parameters of the Weibull distribution given the data in the vector, DATA.  
%
%   [PHAT, PCI] = WEIBFIT(DATA,ALPHA) gives MLEs and 100(1-ALPHA) 
%   percent confidence intervals given the data. By default, the
%   optional parameter ALPHA = 0.05 corresponding to 95% confidence intervals.
%
%   See also WEIBCDF, WEIBINV, WEIBPDF, WEIBPLOT, WEIBRND, WEIBSTAT, MLE. 

%   Reference:
%      [1]  Johnson, Norman L., Kotz, Samuel, & Kemp, Adrienne W.,
%      "Univariate Discrete Distributions, Second Edition", Wiley
%      1992 p. 124-130.

%   Copyright 1993-2011 The MathWorks, Inc.
% $Revision: 1.1.8.2 $  $Date: 2011/07/20 00:09:00 $

warning(message('stats:weibfit:Obsolete'));

if nargin < 2 
    alpha = 0.05;
end
p_int = [alpha/2; 1-alpha/2];

[m,n] = size(data);

if m == 1
   m = n;
   n = 1;
   data = data(:);
end

eprob = [0.5./m:1./m:(m - 0.5)./m]';
w  = log(log(1./(1-eprob)));
z  = sort(log(data));
c  = polyfit(z,w,1);
pinit  = [exp(c(2)) c(1)];
phat = fminsearch('weiblike',pinit,optimset('Display','off'),data);

if nargout == 2
   [logL,info]=weiblike(phat,data);
   sigma = sqrt(diag(info));
   pci = norminv([p_int p_int],[phat; phat],[sigma';sigma']);
end
