function obj = fit(training, group, varargin)
%FIT Create a Naive Bayes classifier object by fitting training data. 
%   NB = NAIVEBAYES.FIT(TRAINING, C) builds a NaiveBayes classifier object
%   NB. TRAINING is an N-by-D numeric matrix of predictor data. Rows of
%   TRAINING correspond to observations; columns correspond to features. C
%   contains the known class labels for TRAINING, and it take one of K
%   possible levels. C is a grouping variable, i.e., it can be a
%   categorical, numeric, or logical vector; a cell vector of strings; or a
%   character matrix with each row representing a class label. Each element
%   of C defines which class the corresponding row of TRAINING belongs to.
%   TRAINING and C must have the same number of rows.
%
%   Type "help groupingvariable" for more information about grouping
%   variables.
%
%   NB = NAIVEBAYES.FIT(..., 'PARAM1',val1, 'PARAM2',val2, ...)
%   specifies one or more of the following name/value pairs:
%
%      'Distribution'  a string or a 1-by-D cell vector of strings,
%                      specifying which distributions FIT uses to model the
%                      data. If the value is a string, FIT models all the
%                      features using one type of distribution. FIT can
%                      also model different features using different types
%                      of distributions. If the value is a cell vector, its
%                      Jth element specifies the distribution FIT uses for
%                      the Jth feature.  The available types of
%                      distributions are:
%          'normal'  (default) Normal (Gaussian) distribution.
%          'kernel'  Kernel smoothing density estimate.
%          'mvmn'    Multivariate multinomial distribution for discrete
%                    data. FIT assumes each individual feature follows a
%                    multinomial model within a class. The parameters for a
%                    feature include the probabilities of all possible
%                    values that the corresponding feature can take.
%          'mn'      Multinomial distribution for classifying the count-
%                    based data such as the bag-of-tokens model. In the
%                    bag-of-tokens model, the value of the Jth feature is
%                    the number of occurrences of the Jth token in this
%                    observation, so it must be a non-negative integer.
%                    When 'mn' is used, FIT considers each observation as
%                    multiple trials of a Multinomial distribution, and
%                    considers each occurrence of a token as one trial.
%                    The number of categories(bins) in this multinomial
%                    model is the number of distinct tokens, i.e., the
%                    number of columns of TRAINING.
%                
%      'Prior'       The prior probabilities for the classes, specified as
%                    one of the following:
%          'empirical'   (default) FIT estimates the prior probabilities
%                        from the relative frequencies of the classes in
%                        TRAINING.
%          'uniform'     The prior probabilities are equal for all classes.
%          vector        A numeric vector of length K specifying the prior
%                        probabilities of the K possible values of C, in
%                        the order described in 'help groupingvariable'.
%          structure     A structure S containing class levels and their
%                        prior probabilities.  S must have two fields:
%                  S.prob  A numeric vector of prior probabilities.
%                  S.group A vector of the same type as C, containing
%                          unique class levels indicating the class for the
%                          corresponding element of prob.
%                        S.group must contain all the K levels in C. It
%                        can also contain classes that do not appear in
%                        C. This can be useful if TRAINING is a subset
%                        of a larger training set. FIT ignores any classes
%                        that appear in S.group but not in C.
%      If the prior probabilities don't sum to one, they will be normalized.
%
%      'KSWidth'     The bandwidth of the kernel smoothing window.  The
%                    default is to select a default bandwidth automatically
%                    for each combination of feature and class, using a
%                    value that is optimal for a Gaussian distribution.
%                    The value can be specified as one of the following:
%          scalar         Width for all features in all classes.
%          row vector     1-by-D vector where the Jth element is the
%                         bandwidth for the Jth feature in all classes.
%          column vector  K-by-1 vector where the Ith element specifies the
%                         bandwidth for all features in the Ith class. K
%                         represents the number of class levels.
%          matrix         K-by-D matrix M where M(I,J) specifies the
%                         bandwidth for the Jth feature in the Ith class.
%          structure      A structure S containing class levels and their
%                         bandwidths.  S must have two fields:
%                  S.width A numeric array of bandwidths specified as a row
%                          vector, or a matrix with D columns.
%                  S.group A vector of the same type as C, containing
%                          unique class levels indicating the class for the
%                          corresponding row of width.
%
%      'KSSupport'   The regions where the density can be applied.  It can
%                    be a string, a two-element vector as shown below, or
%                    a 1-by-D cell array of these values:
%          'unbounded'    (default) The density can extend over the whole
%                         real line.
%          'positive'     The density is restricted to positive values.
%          [L,U]          A two-element vector specifying the finite lower
%                         bound L and upper bound U for the support of the
%                         density.
%
%      'KSType'      The type of kernel smoother to use. It can be a string
%                    or a 1-by-D cell array of strings.  Each string can be
%                    'normal' (default), 'box', 'triangle', or
%                    'epanechnikov'.
%
%  The 'KSWidth', 'KSSupport', and 'KSType' parameters are used only for
%  features with the 'kernel' distribution and are ignored for all others.
%
%  FIT treats NaNs, empty strings or 'undefined' values as missing values.
%  For missing values in C, FIT removes the corresponding rows of
%  TRAINING. For missing values in TRAINING, when distribution 'mn' is
%  used, FIT removes the corresponding rows of TRAINING, otherwise, FIT
%  only removes the missing values and uses the values of other features in
%  the corresponding rows of TRAINING.
%
%  See also NAIVEBAYES, PREDICT, POSTERIOR, PARAMS, GROUPINGVARIABLE,
%  FITDIST, PROBDISTUNIVKERNEL.

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.8.7 $  $Date: 2012/02/02 19:09:45 $

if nargin < 2
    error(message('stats:NaiveBayes:TooFewInputs'));
end

if ~isnumeric(training) 
    error(message('stats:NaiveBayes:fit:TrainingBadType'));
end

if ~isreal(training)
    error(message('stats:NaiveBayes:fit:TrainingComplexType'));
end

obj = NaiveBayes;

[gindex,obj.ClassNames, obj.ClassLevels] = grp2idx(group);
n = size(training,1);
if n == 0 ||size(gindex,1) == 0
    error(message('stats:NaiveBayes:fit:EmptyData'));
end

if n ~= size(gindex,1);
    error(message('stats:NaiveBayes:fit:MismatchedSize'));
end
 
nans = isnan(gindex);
if any(nans)
    training(nans,:) = [];
    gindex(nans) = [];
end

obj.NClasses = length(obj.ClassNames);
obj.ClassSize = hist(gindex,1: obj.NClasses);
obj.CIsNonEmpty = (obj.ClassSize > 0)';
obj.NonEmptyClasses =find(obj.ClassSize>0);

obj.LUsedClasses = length(obj.NonEmptyClasses);
if obj.NClasses > obj.LUsedClasses
    warning(message('stats:NaiveBayes:fit:EmptyGroups'));
end

[n, obj.NDims]= size(training);
if n == 0
    error(message('stats:NaiveBayes:fit:NoData'));
end

% Parse input and error check
pnames = {'distribution' 'prior'   'kswidth'    'kssupport' 'kstype'};
dflts =  {'normal'       'empirical' []         []           []};
[obj.Dist,prior, kernelWidth,obj.KernelSupport, obj.KernelType] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});

if isempty(prior) || (ischar(prior) && strncmpi(prior,'empirical',length(prior)))
    obj.Prior = obj.ClassSize(:)' / sum(obj.ClassSize);
elseif ischar(prior) && strncmpi(prior,'uniform',length(prior))
    obj.Prior = ones(1, obj.NClasses) / obj.NClasses;
    % Explicit prior
elseif isnumeric(prior)
    if ~isvector(prior) || length(prior) ~= obj.NClasses
        error(message('stats:NaiveBayes:fit:ScalarPriorBadSize'));
    end
    obj.Prior = prior;
elseif isstruct(prior)
    if ~isfield(prior,'group') || ~isfield(prior,'prob')
        error(message('stats:NaiveBayes:fit:PriorBadFields'));
    end
    [pgindex,pgroups] = grp2idx(prior.group);
    
    ord =NaN(1,obj.NClasses);
    for i = 1:obj.NClasses
        j = find(strcmp(obj.ClassNames(i), pgroups(pgindex)));
        if isempty(j)
            error(message('stats:NaiveBayes:fit:PriorBadGroup'));
        elseif numel(j) > 1
             error(message('stats:NaiveBayes:fit:PriorDupClasses'));
        else
            ord(i) = j;
        end
    end
    obj.Prior = prior.prob(ord);
    
else
    error(message('stats:NaiveBayes:fit:BadPriorType'));
end

obj.Prior(obj.ClassSize==0) = 0;
if any(obj.Prior < 0) || sum(obj.Prior) == 0
    error(message('stats:NaiveBayes:fit:BadPrior'));
end
obj.Prior = obj.Prior(:)' / sum(obj.Prior); % normalize the row vector

if ischar(obj.Dist)
    obj.Dist = cellstr(obj.Dist); %convert a single string to a cell
elseif ~iscell(obj.Dist)
    error(message('stats:NaiveBayes:fit:BadDist'));
end
% distribution list must be a vector
if ~isvector(obj.Dist)
      error(message('stats:NaiveBayes:fit:BadDist'));
end
if numel(obj.Dist) ~= 1% ~isscalar(obj.Dist)
    if length(obj.Dist) ~= obj.NDims
        error(message('stats:NaiveBayes:fit:BadDistVec'));
    end
    
    try
        u = unique(obj.Dist);
    catch ME
        if isequal(ME.identifier,'MATLAB:UNIQUE:InputClass')
            error(message('stats:NaiveBayes:fit:BadDist'));
        else
            rethrow(ME);
        end
    end
    %if all the distribution are same, make it a scalar cell
    if numel(u) == 1 &&  ~strncmpi(u,'mn',length(u{1}))
        obj.Dist = u;
    end
else
     if ~ischar(obj.Dist{1})
        error(message('stats:NaiveBayes:fit:BadDist'));
     end
end

distNames = {'normal','mvmn','kernel','mn'};
if numel(obj.Dist) == 1 %isscalar(obj.Dist)
    %i = strmatch(lower(obj.Dist),distNames);
    i = find(strncmpi(obj.Dist{1},distNames,length(obj.Dist{1})));
    if isempty(i)
        error(message('stats:NaiveBayes:fit:UnknownScalarDist',obj.Dist{1}));
    elseif numel(i) > 1
        error(message('stats:NaiveBayes:fit:AmbiguousScalarDist',obj.Dist{1}));
    elseif i == 1
        obj.GaussianFS = true(1,obj.NDims);
    elseif i ==2 %'mvmn'
        obj.MVMNFS = true(1,obj.NDims);
    elseif i == 3 %'kernel'
        obj.KernelFS = true(1,obj.NDims);
    end %
    obj.Dist = distNames(i);
    
else %obj.Dist is a vector
    obj.GaussianFS = false(1,obj.NDims); % flag for Gaussian features
    obj.MVMNFS = false(1,obj.NDims); % flag for multivariate multinomial features
    obj.KernelFS = false(1,obj.NDims);   % flag for kernel features

    for d = 1:obj.NDims
        curDist =obj.Dist{d};
        
        i = find(strncmpi(curDist,distNames,length(curDist)));
        if isempty(i)
            error(message('stats:NaiveBayes:fit:UnknownDistVector',curDist,d));
        elseif numel(i)>1
            error(message('stats:NaiveBayes:fit:AmbiguousDistVector',curDist,d));
        elseif i==4
            error(message('stats:NaiveBayes:fit:BadDistMN'));
        elseif i==1
            obj.GaussianFS(d) = true;
        elseif i==2
            obj.MVMNFS(d) = true;
        elseif i==3
            obj.KernelFS(d) = true;
        end
        obj.Dist{d} = distNames{i};
    end %loop over d
    
 u = unique(obj.Dist);
    if length(u) == 1
        obj.Dist = u;
    end
end

if isscalar(obj.Dist) && strcmp(obj.Dist,'mn')
    nans = any(isnan(training),2);%remove rows with any NaN
    %remove rows with invalid values
    trBad =  any(training< 0 |  training ~= round(training), 2);
    if any(trBad)
        warning(message('stats:NaiveBayes:fit:BadDataforMN'));
    end
    t = nans | trBad;
    if any(t)
        training(t,:) = [];
        gindex(t) = [];
    end
else
    nans = all(isnan(training),2);%remove rows with all NaNs
    if any(nans)
        training(nans,:) = [];
        gindex(nans) = [];
    end
    
    for k = obj.NonEmptyClasses
        groupI = (gindex == k);
        if sum(groupI) == 0
            error(message('stats:NaiveBayes:fit:NoDataInEachClass'));
        end
        nanCols =  all(isnan(training(groupI,:)),1);
        if any(nanCols)
            nanCols = strtrim(sprintf('%d ',find(nanCols)));
            error(message('stats:NaiveBayes:fit:TrainingAllNaN', obj.ClassNames{ k }, nanCols));
        end
    end
end

%process the kernel options
if any(obj.KernelFS)
    if ~isempty(kernelWidth)
        if isnumeric(kernelWidth)
            %check the size of kernel width
            [wd1, wd2]=size(kernelWidth);
            if(wd1 ~= 1 && wd1 ~= obj.NClasses) || (wd2 ~= 1 && wd2 ~= obj.NDims)
                error(message('stats:NaiveBayes:fit:ScalarKernelWidthSizeBad'));
            end
            obj.KernelWidth = kernelWidth;
            
        elseif isstruct(kernelWidth)
            if ~isfield(kernelWidth,'group') || ~isfield(kernelWidth,'width')
                error(message('stats:NaiveBayes:fit:KernelWidthBadFields'));
            end
            
            if ~isnumeric(kernelWidth.width)
                error(message('stats:NaiveBayes:fit:KernelWidthNonNumeric'));
            end
            
            [kwgindex,kwgroups] = grp2idx(kernelWidth.group);
            if size(kernelWidth.width,1) ~= length(kwgroups);
                error(message('stats:NaiveBayes:fit:KernelWidthRowSizeBad'));
            end
            if size(kernelWidth.width,2) ~= 1 &&...
                    size(kernelWidth.width,2) ~= obj.NDims;
                error(message('stats:NaiveBayes:fit:KernelWidthColumnSizeBad'));
            end
            ord = NaN(1,obj.NClasses);
            
            for i = 1:obj.NClasses
                j = find(strcmp(obj.ClassNames(i), kwgroups(kwgindex)));
                if isempty(j)
                    error(message('stats:NaiveBayes:fit:KernelWidthBadGroup'));
                elseif numel(j) > 1
                    error(message('stats:NaiveBayes:fit:KernelWidthDupClasses'));
                else
                    ord(i) = j; 
                end
            end
            obj.KernelWidth = kernelWidth.width(ord,:);
        else
            error(message('stats:NaiveBayes:fit:InvalidKernelWidth'));
        end
        
        %check the validity of kernel width.
        if size(obj.KernelWidth,2) > 1
            kwtemp = obj.KernelWidth(:,obj.KernelFS);
        else
            kwtemp = obj.KernelWidth;
        end
        
        if size(obj.KernelWidth,1) > 1
            kwtemp = kwtemp(obj.NonEmptyClasses,:);
        end
        
        kwtemp = kwtemp(:);
        
        if  any(~isfinite(kwtemp)) || any(kwtemp <= 0)
                error(message('stats:NaiveBayes:BadKSWidth'));
        end
        
        
    end % ~isempty(kernelWidth)
    
    if ~isempty(obj.KernelSupport)
        
        if iscell(obj.KernelSupport) 
            if isscalar(obj.KernelSupport) %allow a cell with only one element
                obj.KernelSupport = validSupport(obj.KernelSupport{1});
            else
                if ~isvector(obj.KernelSupport) || length(obj.KernelSupport) ~= obj.NDims
                    error(message('stats:NaiveBayes:fit:BadSupport'));
                end
                %check each kernelsupport
                supporttemp = obj.KernelSupport(obj.KernelFS);
                for i = 1: numel(supporttemp)
                    supporttemp{i}= validSupport(supporttemp{i});
                end
                obj.KernelSupport(obj.KernelFS) = supporttemp;
            end
        else
            obj.KernelSupport = validSupport(obj.KernelSupport);
        end
    else
        obj.KernelSupport = 'unbounded';
    end % ~isempty(obj.KernelSupport)
    
    if ~isempty(obj.KernelType)
        if ischar(obj.KernelType)
            obj.KernelType = cellstr(obj.KernelType);
        elseif ~iscell(obj.KernelType)
            error(message('stats:NaiveBayes:fit:BadKSType'));
        end
        if ~isvector(obj.KernelType)
            error(message('stats:NaiveBayes:fit:BadKSType'));
        end
        
        if isscalar(obj.KernelType)
            obj.KernelType = validKernelType(obj.KernelType{1});
        else
            %check the length of vector kernelType
            if length(obj.KernelType) ~= obj.NDims
                error(message('stats:NaiveBayes:fit:KSTypeBadSize'));
            end
            
            kernelTypeTemp = obj.KernelType(obj.KernelFS);
            for i = 1: numel(kernelTypeTemp)
                kernelTypeTemp{i}= validKernelType(kernelTypeTemp{i});
            end
            obj.KernelType(obj.KernelFS) = kernelTypeTemp;
        end
    else
        obj.KernelType = 'normal';
    end
    
end

obj.Params = cell(obj.NClasses, obj.NDims);

%Start Fit
if isscalar(obj.Dist)
    switch obj.Dist{:}
        case 'mn'
            obj =  mnfit(obj,training, gindex);
        case 'normal'
            obj = gaussianFit(obj, training, gindex);
        case 'mvmn'
            obj = mvmnFit(obj, training,gindex);
        case 'kernel'
            obj = kernelFit(obj,training, gindex);
    end
else
    if any(obj.GaussianFS)
        obj = gaussianFit(obj, training, gindex);
    end
    if any(obj.MVMNFS)
        obj = mvmnFit(obj, training,gindex);
    end
    if any(obj.KernelFS)
        obj = kernelFit(obj,training, gindex);
    end
    
end

end %fit

%--------------------------------------
%estimate parameters using Gaussian distribution
function obj = gaussianFit(obj, training, gidx)
for i = obj.NonEmptyClasses
    groupI = (gidx == i);
    
    gsize = sum(~isnan(training(groupI,obj.GaussianFS)),1);
    if any(gsize < 2)
        error(message('stats:NaiveBayes:fit:NoEnoughDataForGaussian'));
    end
    mu = nanmean(training(groupI,obj.GaussianFS));
    sigma = nanstd(training(groupI,obj.GaussianFS));
    badCols = sigma <= gsize * eps(max(sigma));
    if any(badCols)
        badCols = sprintf('%d ',find(badCols));
        error(message('stats:NaiveBayes:fit:BadVariance', badCols, obj.ClassNames{ i }));
    end
    obj.Params(i,obj.GaussianFS) = mat2cell([mu;sigma],2,...
        ones(1,sum(obj.GaussianFS)));
    %Each cell is a 2-by-1 vector, the first element is the mean,
    %and the second element is the standard deviation.
end
end %function gaussianFit

%-------------------------------------------
%Use kernel density estimate
function obj = kernelFit(obj, training,gidx)

kdfsidx = find(obj.KernelFS);
kw2=[];
if ~isempty(obj.KernelWidth)
    [kwrLen,kwcLen] = size(obj.KernelWidth);
    kw = obj.KernelWidth;
    if kwrLen == 1
        kw = repmat(kw, [obj.NClasses,1]);
    end
    if kwcLen == 1
        kw = repmat(kw, [1,obj.NDims]);
    end
end

for i = obj.NonEmptyClasses
    groupI = (gidx == i);
    
    for j = kdfsidx
        if iscell(obj.KernelSupport)
            kssupport = obj.KernelSupport{j};
        else
            kssupport = obj.KernelSupport;
        end
        if iscell(obj.KernelType)
            kstype = obj.KernelType{j};
        else
            kstype = obj.KernelType;
        end
        
        data = training(groupI,j);
        nans = isnan(data);
        if any(nans)
            data(nans)=[];
            if size(data,1) == 0
                error(message('stats:NaiveBayes:fit:NoEnoughDataForKernel'));
            end
        end
        
        if ~isempty(obj.KernelWidth)
            kw2 = kw(i, j);
        end
        
        obj.Params{i,j} = ...
            fitdist(data,'kernel', 'width',kw2,...
            'support',kssupport,'kernel',kstype);
        
    end
    
    
end
end


%---------------------------
%estimate the parameters using multivariate multinomial
function obj = mvmnFit(obj, training, gidx)

mvmnfsIdx = find(obj.MVMNFS);
d = sum(obj.MVMNFS);
mvmnParams = cell(obj.NClasses,d);
obj.UniqVal = cell(1,d);
for j = 1: d
    data = training(:,mvmnfsIdx(j));
    gidx2 = gidx;
    nans = isnan(data);
    if any(nans)
        data(nans)=[];
        gidx2(nans)=[];
        
    end
    obj.UniqVal{j} = unique(data);
    for i = obj.NonEmptyClasses
        groupI = (gidx2 == i);
        if sum(groupI) == 0
            error(message('stats:NaiveBayes:fit:NoEnoughDataForMVMN'));
        end
        p = histc(data(groupI),obj.UniqVal{j});
        %Add one count for each discrete value of the training data to avoid zero probability
        p= (p+1)/(size(data(groupI),1) +length(obj.UniqVal{1,j}));
        mvmnParams(i,j) = {p(:)};
    end
end

obj.Params(:,obj.MVMNFS) = mvmnParams;

end

%-----------------------------------------------------
% perform Multinomial fit
function obj =  mnfit(obj,training, gidx)
d = size(training,2);
for k = obj.NonEmptyClasses
    groupI = (gidx == k);
    if sum(groupI) == 0
        error(message('stats:NaiveBayes:fit:NoDataInEachClass'));
    end
    
    pw = sum(training(groupI,:),1);
    pw = (pw+1)/(sum(pw)+d);
    %Add one count for  each feature to avoid zero probability
    obj.Params(k,:)= num2cell(pw);% mat2cell(pw,1,ones(1,d));
end
end

%-----------------------------------
%check the validity of kernelsupport
function  kssupport = validSupport(kssupport)

badSupport = false;
if ischar(kssupport) && size(kssupport,1)==1
    supportName = {'unbounded' 'positive'};
  
    i = find(strncmpi(kssupport,supportName,length(kssupport)));
    if isempty(i)
        badSupport = true;
    else
        kssupport = supportName{i};
    end
    
elseif ~(isnumeric(kssupport) && numel(kssupport)==2 ...
        && all(isfinite(kssupport)) && kssupport(1) < kssupport(2))
    badSupport = true;
end
if badSupport
    error(message('stats:NaiveBayes:fit:BadSupport'));
end
end

%----------------------------------------------------------------
%check the validity of kernel Type
function type = validKernelType(type)
typeNames ={'normal' , 'box', 'triangle', 'epanechnikov'};

if ~ischar(type)
    error(message('stats:NaiveBayes:fit:BadKSType'));
end

i = find(strncmpi(type,typeNames, length(type)));
if isempty(i)
    error(message('stats:NaiveBayes:fit:UnknownKSType',type));   
end
type= typeNames{i};
end
