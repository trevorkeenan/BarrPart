function name = dfgetfitname()

% $Revision: 1.1.8.2 $  $Date: 2011/08/17 21:29:43 $
% Copyright 2003-2011 The MathWorks, Inc.

count=dfgetset('fitcount');
if isempty(count)
    count = 1;
end
taken = 1;
while taken
    name=getString(message('stats:dfstrings:name_FitNumber', count));
    if isempty(find(getfitdb,'name',name))
        taken = 0;
    else
        count=count+1;
    end
end
dfgetset('fitcount',count+1);

