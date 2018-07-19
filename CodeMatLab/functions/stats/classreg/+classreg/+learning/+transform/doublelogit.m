function out = doublelogit(in)

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2010/09/20 15:09:41 $

out = 1./(1+exp(-2*in));
end
