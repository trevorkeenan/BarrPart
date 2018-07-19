function out = symmetricismax(in)

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2010/09/20 15:09:50 $

out = classreg.learning.transform.ismax(in);

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2010/09/20 15:09:50 $
out(out<0.5) = -1;
end
