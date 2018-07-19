function [postp,cidx,logP] = posterior(obj,test,varargin)
%POSTERIOR Compute posterior probability of each class for test data.
%    POST = POSTERIOR(NB,TEST) returns the posterior probability of the
%    observations in TEST according to the NaiveBayes object NB. TEST is a
%    N-by-NB.NDIMS matrix, where N is the number of observations in the
%    test data. Rows of TEST correspond to points, columns of TEST
%    correspond to features. POST is a N-by-NB.NCLASSES matrix containing
%    the posterior probability of each observation for each class.
%    POST(I,J) is the posterior probability of point I belonging to class
%    J. Classes are ordered the same as NB.CLASSLEVELS, i.e., column J of
%    POST corresponds to the Jth class in NB.CLASSLEVELS. The posterior
%    probabilities corresponding to any empty classes are NaN.
%
%    [POST,CPRE] = POSTERIOR(NB,TEST) returns CPRE, an N-by-1 vector,
%    containing the class to which each row of TEST has been assigned.
%    CPRE has the same type as NB.CLASSLEVELS.
%
%    [POST,CPRE,LOGP] = POSTERIOR(NB,TEST) returns LOGP, an N-by-1
%    vector containing estimates of the log of the probability density
%    function (PDF). LOGP(I) is the log of the PDF of point I. The PDF
%    value of point I is the sum of
%
%         Prob(point I | class J) * Pr{class J}
%
%    taken over all classes.
%
%    [...] = POSTERIOR(..., 'HandleMissing',VAL) specifies how NaN (missing
%    values) are treated. VAL can be one of the following:
%       'off' (default) Observations with NaN in any columns are not
%             classified into any class. The corresponding rows in POST and
%             LOGP are NaN. The corresponding rows in CPRE are NaN (if
%             OBJ.CLASSLEVELS is numeric or logical), empty strings (if
%             OBJ.CLASSLEVELS is char or cell array of strings) or <undefined>
%             (if OBJ.CLASSLEVELS is categorical).
%       'on'  For observations having NaN in some (but not all) columns,
%             POST and CPRE are computed using the columns with non-NaN
%             values.  Corresponding LOGP values are NaN.
%
%   See also NAIVEBAYES, FIT, PREDICT.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2011/05/09 01:27:41 $

if nargin < 2
    error(message('stats:NaiveBayes:posterior:TooFewInputs'));
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

wasAnyNans = any(isnan(test),2);
if strcmp(handleMissing,'off')
    wasInvalid = wasAnyNans;
else
    wasInvalid = all(isnan(test),2);
end

if isscalar(obj.Dist) && strcmp(obj.Dist,'mn')
    testBad = any(test<0 |test ~= round(test),2);
    if any(testBad)
         warning(message('stats:NaiveBayes:posterior:BadDataforMN'));
       wasInvalid = wasInvalid | testBad;
    end
end

hadInvalid = any(wasInvalid);
if hadInvalid
    test(wasInvalid,:)= [];
end

log_condPdf = getlogCondPDF(obj, test, handleMissing);

if nargout >= 3
    [cidx, postp, logP] = getClassIdx(obj,log_condPdf);
    %convert  class index to the corresponding class levels
    cidx = cidx2glevel(cidx, obj.ClassLevels);
    
    if hadInvalid
        [cidx,postp,logP] = dfswitchyard('statinsertnan',wasInvalid,cidx,postp,logP);
    end
    if strcmp(handleMissing,'on')
        logP(wasAnyNans) = NaN;
    end
else
    [cidx,postp] = getClassIdx(obj,log_condPdf);
    %convert integer class index to the corresponding class levels
    cidx = cidx2glevel(cidx,obj.ClassLevels);
    
    if hadInvalid
        [cidx,postp] = dfswitchyard('statinsertnan',wasInvalid,cidx,postp);
    end
end

end %posterior

%--------------------
%handle NaN values in cidx
function cidx = cidx2glevel(cidx,ClassLevels)

isCidxNaN= isnan(cidx);
cidx = ClassLevels(cidx(~isCidxNaN),:);
if any(isCidxNaN)
    cidx= dfswitchyard('statinsertnan',isCidxNaN,cidx);
end
end

