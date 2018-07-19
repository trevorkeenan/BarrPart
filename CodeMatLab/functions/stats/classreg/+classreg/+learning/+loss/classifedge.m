function e = classifedge(C,Sfit,W,cost)

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2010/09/20 15:09:08 $


m = classreg.learning.loss.classifmargin(C,Sfit);

% Average the margins
notNaN = ~isnan(m);
e = sum(W(notNaN).*m(notNaN)) / sum(W(notNaN));
end
