function out = ismax(in)

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2010/09/20 15:09:45 $

[N,K] = size(in);
out = zeros(N,K);
[~,colnum] = max(in,[],2);
for k=1:K
    out(colnum==k,k) = 1;
end
end
