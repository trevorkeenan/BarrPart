classdef FullClassificationRegressionModel < classreg.learning.Predictor
%FullClassificationRegressionModel Full classification or regression model.
%   FullClassificationRegressionModel is the super class for full
%   classification or regression models represented by objects storing the
%   training data.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.8 $  $Date: 2012/01/27 21:05:21 $

    properties(GetAccess=public,SetAccess=protected)
        %X X data used to train this model.
        %   The X property is a numeric matrix of size N-by-P, where N is the
        %   number of observations (rows) and P is the number of predictors
        %   (columns) in the training data.  This matrix contains the predictor
        %   values.
        %
        %   See also classreg.learning.FullClassificationRegressionModel.
        X = [];
        
        %W Weights of observations used to train this model.
        %   The W property is a numeric vector of size N, where N is the
        %   number of observations. The sum of weights is 1.
        %
        %   See also classreg.learning.FullClassificationRegressionModel.
        W = [];
        
        %MODELPARAMS Model parameters.
        %   The ModelParams property holds parameters used for training this model.
        %
        %   See also classreg.learning.FullClassificationRegressionModel.
        ModelParams = [];
    end

    properties(GetAccess=public,SetAccess=protected,Hidden=true)
        PrivY = [];
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %NOBSERVATIONS Number of observations.
        %   The NObservations property is a numeric scalar holding the number of
        %   observations in the training data.
        %
        %   See also classreg.learning.FullClassificationRegressionModel.
        NObservations;
    end
    
    methods
        function n = get.NObservations(this)
            n = size(this.X,1);
        end
    end

    methods(Abstract,Static)
        obj = fit(X,Y,varargin)
    end
    
    methods(Abstract)
        cmp = compact(this)
        partModel = crossval(this,varargin)
    end
    
    methods(Access=protected)
        function this = FullClassificationRegressionModel(dataSummary,X,Y,W,modelParams)
            this = this@classreg.learning.Predictor(dataSummary);
            this.X = X;
            this.PrivY = Y;
            this.W = W;
            this.ModelParams = modelParams;
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.Predictor(this,s);
            s.NObservations = this.NObservations;
        end
    end
    
    methods(Static,Hidden)
        function [X,Y,W,dataSummary] = prepareDataCR(X,Y,varargin)
            % Process input args
            args = {'weights' 'predictornames' 'responsename' 'categoricalpredictors'};
            defs = {       []               []             []                      []};
            [W,predictornames,responsename,catpreds] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check input type
            if ~isnumeric(X) || ndims(X)>2
                error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:BadXType'));
            end
            
            % Check input size
            if isempty(X) || isempty(Y)
                error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:NoData'));
            end
            N = size(X,1);
            if N~=length(Y)
                error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:InputSizeMismatch'));
            end
 
            % Check weights and normalize to 1
            if isempty(W)
                W = ones(N,1);
            else
                if ~isfloat(W) || length(W)~=size(X,1) || ~isvector(W)
                    error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:BadW'));
                end
                if any(W<0) || all(W==0)
                    error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:NegativeWeights'));
                end
                W = W(:);
            end

            % Get rid of instances that have NaN's in all predictors
            t = all(isnan(X),2);
            if any(t)
                Y(t) = [];
                X(t,:) = [];
                W(t) = [];
            end
            if isempty(X)
                error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:NoGoodXData'));
            end
            
            % Get rid of observations with zero weights or NaNs
            t = (W==0 | isnan(W));
            if any(t)
                Y(t) = [];
                X(t,:) = [];
                W(t) = [];
            end
            if isempty(X)
                error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:NoGoodWeights'));
            end         

            % Process predictor names
            D = size(X,2);
            if     isempty(predictornames)
                predictornames = D;
            elseif isnumeric(predictornames)
                if ~(isscalar(predictornames) && predictornames==D)
                    error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:BadNumericPredictor', D));
                end
            else
                if ~iscellstr(predictornames)
                    if ~ischar(predictornames)
                        error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:BadPredictorType'));
                    end
                    predictornames = cellstr(predictornames);
                end
                if length(predictornames)~=D
                    error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:PredictorMismatch', D));
                end
            end
            predictornames = predictornames(:)';
            
            % Process response name
            if isempty(responsename)
                responsename = 'Y';
            else
                if ~ischar(responsename)
                    error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:BadResponseName'));
                end
            end
                
            % Find categorical predictors
            if isnumeric(catpreds) % indices of categorical predictors
                catpreds = ceil(catpreds);
                if any(catpreds<1) || any(catpreds>D)
                    error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:BadCatPredIntegerIndex', D));
                end
            elseif islogical(catpreds)
                if length(catpreds)~=D
                    error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:BadCatPredLogicalIndex', D));
                end
                idx = 1:D;
                catpreds = idx(catpreds);
            elseif ischar(catpreds) && strcmpi(catpreds,'all')
                catpreds = 1:D;
            else
                if ~ischar(catpreds) && ~iscellstr(catpreds)
                    error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:BadCatVarType'));
                end
                if ~iscellstr(catpreds)
                    catpreds = cellstr(catpreds);
                end
                if isnumeric(predictornames)
                    error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:CharCatVarWithoutVarNames'));
                end
                [tf,pos] = ismember(catpreds,predictornames);
                if any(~tf)
                    error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:BadCatVarName', catpreds{ find( tf, 1, 'first' ) }));
                end
                catpreds = pos;
            end
            
            % Summarize data
            dataSummary.PredictorNames = predictornames;
            dataSummary.CategoricalPredictors = catpreds;
            dataSummary.ResponseName = responsename;
        end
        
        function catchWeights(varargin)
            args = {'weights'};
            defs = {       []};
            [w,~,~] = internal.stats.parseArgs(args,defs,varargin{:});
            if ~isempty(w)
                error(message('stats:classreg:learning:FullClassificationRegressionModel:catchWeights:NonEmptyWeights'));
            end
        end

    end
    
end
