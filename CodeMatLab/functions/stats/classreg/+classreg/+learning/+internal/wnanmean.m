function M = wnanmean(X,W)
%WNANMEAN Weighted mean values ignoring NaN's.
%   M=WNANMEAN(X,W) returns mean for matrix X along 1st dimension weighted
%   by vector W.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2011/03/22 18:38:28 $

W(isnan(W)) = 0;
W = W(:);

X = bsxfun(@times,X,W);
tfnan = isnan(X);
X(tfnan) = 0;

Wcol = sum(bsxfun(@times,~tfnan,W),1);
M = sum(X,1) ./ Wcol;
end
