function e = classiferror(C,Sfit,W,~)

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2011/11/09 17:46:16 $

% Classification error is the fraction of incorrect predictions
notNaN = ~all(isnan(Sfit),2);
[~,y] = max(C(notNaN,:),[],2);
[~,yfit] = max(Sfit(notNaN,:),[],2);
W = W(notNaN);
e = sum((y~=yfit).*W) / sum(W);
end
