function out = sign(in)

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2010/09/20 15:09:47 $

out = zeros(size(in));
out(in<0) = -1;
out(in>0) = +1;
end
