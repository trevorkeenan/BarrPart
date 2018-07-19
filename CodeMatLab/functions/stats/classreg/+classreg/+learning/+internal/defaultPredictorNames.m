function names = defaultPredictorNames(D)

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2011/11/09 17:46:14 $

names = strcat({'x'},num2str((1:D)','%-d'))';
end
