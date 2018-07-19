function m = classifmargin(C,Sfit)

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2010/09/20 15:09:10 $

[N K] = size(C);

if K==1
    m = NaN(N,1);
    return;
end

[~,trueC] = max(C,[],2);
m = zeros(N,1);
for k=1:K
    trueK = false(K,1);
    trueK(k) = true;
    idx = trueC==k;
    m(idx) = Sfit(idx,trueK) - max(Sfit(idx,~trueK),[],2);
end
end
