function loss = mse(Y,Yfit,W)

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2010/09/20 15:09:12 $

notNaN = ~isnan(Yfit);
loss = sum(W(notNaN) .* (Y(notNaN)-Yfit(notNaN)).^2) / sum(W(notNaN));
end
