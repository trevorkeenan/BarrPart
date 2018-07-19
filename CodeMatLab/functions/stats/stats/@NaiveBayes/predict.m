function gpre = predict(obj,test,varargin)
%PREDICT Predict class label for test data.
%   CPRE = PREDICT(NB,TEST) classifies each row in TEST into one of
%   the classes according to the NaiveBayes classifier NB, and returns the
%   predicted class level CPRE. TEST is an N-by-NB.NDIMS matrix, where N
%   is the number of observations in the test data. Rows of TEST correspond
%   to points, columns of TEST correspond to features. CPRE is an N-by-1
%   vector of the same type as NB.CLASSLEVELS, and it indicates the class
%   to which each row of TEST has been assigned.
%
%   CPRE = PREDICT(..., 'HandleMissing',VAL) specifies how NaN (missing
%   values) in TEST are treated. VAL can be one of the following:
%       'off' (default) Observations with NaN in any of columns are not
%             classified into any class. The corresponding rows in CPRE are
%             NaN (if OBJ.CLASSLEVELS is numeric or logical), empty strings
%             (if OBJ.CLASSLEVELS is char or cell array of strings) or
%             <undefined> (if OBJ.CLASSLEVELS is categorical).
%       'on'  For observations having NaN values in some (but not all) of
%             columns CPRE is computed using the columns with non-NaN
%             values. 
%
%   See also NAIVEBAYES, FIT, POSTERIOR.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2011/05/09 01:27:42 $

if nargin < 2
    error(message('stats:NaiveBayes:predict:TooFewInputs'));
end

if ~isnumeric(test) 
    error(message('stats:NaiveBayes:predict:TestBadType'));
end

if ~isreal(test)
    error(message('stats:NaiveBayes:predict:TestComplexType'));
end

if size(test,2)~= obj.NDims
      error(message('stats:NaiveBayes:BadTestSize', obj.NDims));
end
 
pnames = {'handlemissing'};
dflts = {'off'};
handleMissing = internal.stats.parseArgs(pnames, dflts, varargin{:});
handleMissing = internal.stats.getParamVal(handleMissing,{'on' 'off'},'HandleMissing');

if strcmp(handleMissing,'off')
    wasInvalid =  any(isnan(test),2);
else
    wasInvalid = all(isnan(test),2);
end

if isscalar(obj.Dist) && strcmp(obj.Dist,'mn')
    testBad = any(test<0 |test ~= round(test),2);
    if any(testBad)
        warning(message('stats:NaiveBayes:predict:BadDataforMN'));
         wasInvalid = wasInvalid | testBad;
    end
end

hadInvalid = any(wasInvalid);
if hadInvalid
    test(wasInvalid,:)= [];
end

log_condPdf = getlogCondPDF(obj,test,handleMissing);
gpre = getClassIdx(obj,log_condPdf);
%convert class index to the corresponding class levels

isGpreNaN= isnan(gpre);
gpre= obj.ClassLevels(gpre(~isGpreNaN),:);
if any(isGpreNaN)
    try
      gpre= dfswitchyard('statinsertnan',isGpreNaN,gpre);
    catch ME
        if ~isequal(ME.identifier,'stats:statinsertnan:InputTypeIncorrect')
                 rethrow(ME);
             else
                 error(message('stats:NaiveBayes:predict:logicalwithNaN'));
         end
    end
end

if hadInvalid
    try 
    gpre = dfswitchyard('statinsertnan', wasInvalid, gpre);
    catch ME
           if ~isequal(ME.identifier,'stats:statinsertnan:InputTypeIncorrect')
                 rethrow(ME);
             else
                 error(message('stats:NaiveBayes:predict:logicalwithNaN'));
           end
    end
end

end
