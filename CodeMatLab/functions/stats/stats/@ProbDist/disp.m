function disp(obj)
%DISP Display a probability distribution object.
%   DISP(PD) prints a text representation of the probability distribution
%   PD, without printing the object name.  In all other ways it's
%   the same as leaving the semicolon off an expression.
%
%   See also ProbDist, ProbDist/DISPLAY.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2011/08/17 21:29:18 $

isLoose = strcmp(get(0,'FormatSpacing'),'loose');

if (isLoose)
    fprintf('\n');
end

fprintf(getString(message('stats:ProbDist:disp_Distribution',obj.DistName)));

if (isLoose)
    fprintf('\n');
end
