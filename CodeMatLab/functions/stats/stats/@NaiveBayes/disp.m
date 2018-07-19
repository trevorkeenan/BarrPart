function disp(obj)
%DISP Display a NAIVEBAYES classifier object.
%   DISP(NB) prints a text representation of the NaiveBayes object NB, without
%   printing the object name.  In all other ways it's the same as leaving
%   the semicolon off an expression.
%
%   See also NAIVEBAYES, NAIVEBAYES/DISPLAY.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2011/10/15 02:07:23 $

isLoose = strcmp(get(0,'FormatSpacing'),'loose');
if (isLoose)
    fprintf('\n');
end

fprintf('%s\n',getString(message('stats:NaiveBayes:disp_NBClassifier',...
    obj.LUsedClasses,obj.NDims)));
if length(obj.Dist) < 10
   fprintf(getString(message('stats:NaiveBayes:disp_FeatureDistributions')))
    for i = 1:length(obj.Dist)-1
       fprintf('%s, ', obj.Dist{i});
    end
   fprintf('%s\n',obj.Dist{end});
end

if obj.LUsedClasses < 10
    fprintf(getString(message('stats:NaiveBayes:disp_Classes')))
    for i = obj.NonEmptyClasses(1:end-1)
       fprintf('%s, ', obj.ClassNames{i});
    end
    fprintf('%s\n', obj.ClassNames{obj.NonEmptyClasses(end)});
end

if (isLoose)
    fprintf('\n');
end
