function r = gammaincratio(x,s)
% GAMMAINCRATIO Ratio of incomplete gamma function values at S and S-1.
%    R=GAMMAINCRATIO(X,S) computes GAMMAINC(X,S)/GAMMAINC(X,S-1). S
%    must be greater or equal to 1. X and S must be of the same size.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2010/09/20 15:10:03 $

% Initialize
r = zeros(size(s));

% What is small S?
smalls = s<2 | s<=x;
if any(~smalls(:))
    idx = find(~smalls);
    smalls(idx) = s(idx).*log(x(idx))-gammaln(s(idx)+1) > log(realmin(class(x)));
end

% For small S, use the ratio computed directly
if any(smalls(:))
    r(smalls) = gammainc(x(smalls),s(smalls))./gammainc(x(smalls),s(smalls)-1);
end

% For large S, estimate numerator and denominator using QUADL.
if any(~smalls(:))
    idx = find(~smalls);
    x = x(idx);
    s = s(idx);
    for i=1:numel(x)
        f1 = @(t) s(i)*t.^(s(i)-1).*exp(x(i)*(1-t));
        f0 = @(t) (s(i)-1)*t.^(s(i)-2).*exp(x(i)*(1-t));
        r(idx(i)) = x(i)/s(i)*quadl(f1,0,1)/quadl(f0,0,1);
    end
end
end
