function funloss = lossCheck(funloss,type)  

%   Copyright 2010-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2011/11/09 17:46:15 $  

if     ischar(funloss)
    if     strcmp(type,'classification')
        allowed = {'binodeviance' 'classifedge' 'classiferror' 'classifmargin' 'exponential' 'mincost'};
    elseif strcmp(type,'regression')
        allowed = {'mse'};
    else
        allowed = {};
    end
    idx = find(strncmpi(funloss,allowed,length(funloss)));
    if isempty(idx) || ~isscalar(idx)
        error(message('stats:classreg:learning:internal:lossCheck:BadFunlossString'));
    end
    funloss = str2func(['classreg.learning.loss.' allowed{idx}]);
elseif ~isa(funloss,'function_handle')
    error(message('stats:classreg:learning:internal:lossCheck:BadFunlossType'));
end
end
