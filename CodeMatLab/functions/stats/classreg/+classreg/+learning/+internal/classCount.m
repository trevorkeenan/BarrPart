function C = classCount(expectedY,observedY)

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2011/05/09 01:23:40 $

if isempty(observedY)
    C = false(0);
    return;
end
[tf,grp] = ismember(observedY,expectedY);
if ~all(tf)
    idx = find(~tf,1,'first');
    if     isa(observedY,'classreg.learning.internal.ClassLabel') ...
            || isa(observedY,'categorical') || iscellstr(observedY)
        str = char(observedY(idx));
    else
        str = num2str(observedY(idx,:));
    end
    if isa(observedY,'classreg.learning.internal.ClassLabel')
        cls = class(labels(observedY));
    else
        cls = class(observedY);
    end
    error(message('stats:classreg:learning:internal:classCount:UnknownClass', str, cls));
end
N = length(observedY);
K = length(expectedY);
C = false(N,K);
C(sub2ind([N K],(1:N)',grp)) = true;
end
