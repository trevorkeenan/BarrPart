function loss = exponential(C,Sfit,W,cost)

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2010/09/20 15:09:11 $


[~,y] = max(C,[],2);
[N,K] = size(C);
Sfit = Sfit(sub2ind([N K],(1:N)',y));
notNaN = ~isnan(Sfit);
loss = sum( W(notNaN).*exp(-Sfit(notNaN)) ) / sum(W(notNaN));
end
