function model = stepwise(X,varargin) % [X, y | DS], start, ...
%STEPWISE Create linear regression model by stepwise regression.
%   LM = LinearModel.stepwise(D,MODELSPEC) fits the model specified by
%   MODELSPEC to variables in the dataset D, adds or removes terms by
%   stepwise regression, and returns the linear model LM. MODELSPEC can be
%   any of the values accepted by the LinearModel.fit function. The default
%   is 'constant' to start with only the constant term in the model.
%
%   After the initial fit, the function examines a set of available terms,
%   and adds the best one to the model if an F-test for adding the term has
%   a p-value 0.05 or less. If no terms can be added, it examines the terms
%   currently in the model, and removes the worst one if an F-test for
%   removing it has a p-value 0.10 or greater. It repeats this process
%   until no terms can be added or removed. The function never removes the
%   constant term. It may add terms defined by MODELSPEC, as well as terms
%   that are two-way interactions among the predictors.
%
%   LM = LinearModel.stepwise(X,Y) fits a linear regression model using the
%   column vector Y as a response variable and the columns of the matrix X
%   as predictor variables, performs stepwise regression, and returns the
%   final result as the linear model LM.
%
%   LM = LinearModel.stepwise(X,Y,MODELSPEC) uses the model specified by
%   MODELSPEC as the initial model. See LinearModel.fit for valid MODELSPEC
%   values.
%
%   LM = LinearModel.stepwise(...,PARAM1,VAL1,PARAM2,VAL2,...) specifies
%   one or more of the following name/value pairs:
%
%      'Weights'          Vector of N non-negative weights, where N is the
%                         number of rows in DS or Y. Default is ones(N,1).
%      'VarNames'         Cell array of strings specifying the names to use
%                         for the columns in X. Default is {'x1','x2',...}
%                         for the predictors and 'y' for the response.
%                         Not allowed when fitting to a dataset.
%      'CategoricalVars'  Vector of integer or logical indices specifying
%                         the variables in DS or the columns in X that
%                         should be treated as categorical. Default is to
%                         treat DS variables as categorical if they are
%                         categorical, logical, or char arrays, or cell
%                         arrays of strings.
%      'Exclude'          Vector of integer or logical indices into the
%                         rows of DS, or X and Y, that should be excluded
%                         from the fit. Default is to use all rows.
%      'Intercept'        true (default) to include a constant term in the
%                         model, or false to omit it.
%      'PredictorVars'    A specification of the variables to use as
%                         predictors, either as a cell array of variable
%                         names, or a vector of integer or logical indices
%                         into the variables in DS or the columns in X.
%      'ResponseVar'      The response variable, specified either as a
%                         variable name or number.
%      'Lower'            A model specification defining the terms that
%                         cannot be removed from the model. Default
%                         'constant', meaning only the intercept.
%      'Upper'            A model specification defining the terms
%                         available to be added to the model. Default
%                         'interactions' for pairwise interaction terms.
%      'Criterion'        The criterion to be used in choosing terms to add
%                         or remove, chosen from 'SSE' (default), 'AIC',
%                         'BIC', 'Rsquared', 'AdjRsquared'.
%      'PEnter'           For the 'SSE' criterion, a value E such that a
%                         term may be added if its p-value is less than or
%                         equal to E. For the other criteria, a term may be
%                         added if the improvement in the criterion is at
%                         least E.
%      'PRemove'          For the 'SSE' criterion, a value R such that a
%                         term may be removed if its p-value is greater or
%                         equal to R. For the other criteria, a term may be
%                         added if it reduces the criterion no more than R.
%      'NSteps'           The maximum number of steps that may be taken,
%                         starting from the initial model. Default is no
%                         limit on the number of steps.
%      'Verbose'          An integer from 0 to 2 controlling the display of
%                         information. Verbose=1 (the default) displays the
%                         action taken at each step. Verbose=2 also
%                         displays the actions evaluated at each step.
%                         Verbose=0 suppresses all display.
%
%   The following table shows the default 'PEnter' and 'PRemove' values for
%   the different criteria, and indicates which must be larger than the
%   other:
%
%      Criterion     PEnter   PRemove    Compared against
%      'SSE'         0.05   < 0.10       p-value for F test
%      'AIC'         0      < 0.01       change in AIC
%      'BIC'         0      < 0.01       change in BIC
%      'Rsquared'    0.1    > 0.05       increase in R-squared
%      'AdjRsquared' 0      > -0.05      increase in adjusted R-squared
%
%   Example:
%      % Perform stepwise regression, starting from the constant model
%      % (intercept only), adding linear terms as required.
%      load hald
%      lm = LinearModel.stepwise(ingredients,heat,'constant','upper','linear')
%
%   See also LinearModel, LinearModel.fit.

%   Copyright 2011-2012 The MathWorks, Inc.

[X,y,haveDataset,otherArgs] = LinearModel.handleDataArgs(X,varargin{:});


paramNames = {'Intercept' 'PredictorVars' 'ResponseVar' 'Weights' 'Exclude' 'CategoricalVars' ...
    'VarNames' 'Lower' 'Upper' 'Criterion' 'PEnter' 'PRemove' 'NSteps' 'Verbose'};
paramDflts = {true [] [] [] [] [] [] 'constant' 'interactions' 'SSE' [] [] Inf 1};

% Default model is constant only.
if isempty(otherArgs)
    start = 'constant';
else
    arg1 = otherArgs{1};
    if mod(length(otherArgs),2)==1 % odd, model followed by pairs
        start = arg1;
        otherArgs(1) = [];
    elseif internal.stats.isString(arg1) && ...
            any(strncmpi(arg1,paramNames,length(arg1)))
        % omitted model but included name/value pairs
        start = 'constant';
    end
end

[intercept,predictorVars,responseVar,weights,exclude,asCatVar, ...
    varNames,lower,upper,crit,penter,premove,nsteps,verbose,supplied] = ...
    internal.stats.parseArgs(paramNames, paramDflts, otherArgs{:});

[penter,premove] = classreg.regr.TermsRegression.getDefaultThresholds(crit,penter,premove);

if ~isscalar(verbose) || ~ismember(verbose,0:2)
    error(message('stats:LinearModel:BadVerbose'));
end

% *** need more reconciliation among start, lower, upper, and between those
% and intercept, predictorVars, varNames

if ~supplied.ResponseVar && (classreg.regr.LinearFormula.isTermsMatrix(start) || classreg.regr.LinearFormula.isModelAlias(start))
    if isa(lower,'classreg.regr.LinearFormula')
        responseVar = lower.ResponseName;
        supplied.ResponseVar = true;
    else
        if internal.stats.isString(lower) && ~classreg.regr.LinearFormula.isModelAlias(lower)
            lower = LinearModel.createFormula(supplied,lower,X, ...
                        predictorVars,responseVar,intercept,varNames,haveDataset);
            responseVar = lower.ResponseName;
            supplied.ResponseVar = true;
        elseif isa(upper,'classreg.regr.LinearFormula')
            responseVar = upper.ResponseName;
            supplied.ResponseVar = true;
        else
            if internal.stats.isString(upper) && ~classreg.regr.LinearFormula.isModelAlias(upper)
                upper = LinearModel.createFormula(supplied,upper,X, ...
                            predictorVars,responseVar,intercept,varNames,haveDataset);
                responseVar = upper.ResponseName;
                supplied.ResponseVar = true;
            end
        end
    end
end
    
if ~isa(start,'classreg.regr.LinearFormula')
    ismodelalias = classreg.regr.LinearFormula.isModelAlias(start);
    start = LinearModel.createFormula(supplied,start,X, ...
                predictorVars,responseVar,intercept,varNames,haveDataset);
else
    ismodelalias = false;
end

if ~isa(lower,'classreg.regr.LinearFormula')
    if classreg.regr.LinearFormula.isModelAlias(lower)
        if supplied.PredictorVars
            lower = {lower,predictorVars};
        end
    end
    lower = classreg.regr.LinearFormula(lower,start.VariableNames,start.ResponseName,start.HasIntercept,start.Link);
end
if ~isa(upper,'classreg.regr.LinearFormula')
    if classreg.regr.LinearFormula.isModelAlias(upper)
        if supplied.PredictorVars
            upper = {upper,predictorVars};
        end
    end
    upper = classreg.regr.LinearFormula(upper,start.VariableNames,start.ResponseName,start.HasIntercept,start.Link);
end

% Remove categorical powers, if any, but warning only if these powers were
% requested explicitly, not just created via something like 'quadratic'
nvars = size(X,2);
if haveDataset
    isCat = datasetfun(@internal.stats.isDiscreteVar,X);
else
    isCat = [repmat(internal.stats.isDiscreteVar(X),1,nvars) internal.stats.isDiscreteVar(y)];
    nvars = nvars+1;
end
if ~isempty(asCatVar)
    isCat = classreg.regr.FitObject.checkAsCat(isCat,asCatVar,nvars,haveDataset,start.VariableNames);
end
if any(isCat)
    start = removeCategoricalPowers(start,isCat,ismodelalias);
    lower = removeCategoricalPowers(lower,isCat,ismodelalias);
    upper = removeCategoricalPowers(upper,isCat,ismodelalias);
end

if haveDataset
    model = LinearModel.fit(X,start.Terms,'ResponseVar',start.ResponseName, ...
        'Weights',weights,'Exclude',exclude,'CategoricalVars',asCatVar,'RankWarn',false);
else
    model = LinearModel.fit(X,y,start.Terms,'ResponseVar',start.ResponseName, ...
        'Weights',weights,'Exclude',exclude,'CategoricalVars',asCatVar, ...
        'VarNames',start.VariableNames,'RankWarn',false);
end

model.Steps.Start = start;
model.Steps.Lower = lower;
model.Steps.Upper = upper;
model.Steps.Criterion = crit;
model.Steps.PEnter = penter;
model.Steps.PRemove = premove;
model.Steps.History = [];

model = stepwiseFitter(model,nsteps,verbose);
checkDesignRank(model);
end
