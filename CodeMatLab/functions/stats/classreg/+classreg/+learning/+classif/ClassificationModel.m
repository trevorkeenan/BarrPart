classdef ClassificationModel < classreg.learning.Predictor
%ClassificationModel Compact classification model.
%   ClassificationModel is the super class for compact classification
%   models.

%   Copyright 2010-2012 The MathWorks, Inc.
%   $Revision: 1.1.6.13 $  $Date: 2012/01/27 21:05:25 $

    properties(GetAccess=public,SetAccess=protected,Hidden=true)
        ClassSummary = struct('ClassNames',{},'NonzeroProbClasses',{},'Cost',[],'Prior',[]);
    end
    
    properties(GetAccess=public,SetAccess=protected,Hidden=true)
        PrivScoreTransform = [];
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %CLASSNAMES Names of classes in Y.
        %   The ClassNames property is an array containing the class names for the
        %   response variable Y.
        %
        %   See also classreg.learning.classif.ClassificationModel.
        ClassNames;
    end

    properties(GetAccess=public,SetAccess=public,Dependent=true)
        %PRIOR Prior class probabilities.
        %   The Prior property is a vector with prior probabilities for classes.
        %
        %   See also classreg.learning.classif.ClassificationModel.
        Prior;
        
        %COST Misclassification costs.
        %   The Cost property is a square matrix with misclassification costs.
        %   Cost(I,J) is the cost of misclassifying class ClassNames(I) as class
        %   ClassNames(J).
        %
        %   See also classreg.learning.classif.ClassificationModel.
        Cost;

        %SCORETRANSFORM Transformation applied to predicted classification scores.
        %   The ScoreTransform property is a string describing how raw
        %   classification scores predicted by the model are transformed. You can
        %   assign a function handle or one of the following strings to this
        %   property: 'none', 'doublelogit', 'identity', 'invlogit', 'ismax',
        %   'logit', 'sign', 'symmetricismax', 'symmetriclogit', and 'symmetric'.
        %   You can use either 'identity' or 'none' for the identity
        %   transformation.
        %
        %   See also classreg.learning.classif.ClassificationModel.
        ScoreTransform;
    end
    
    properties(GetAccess=public,SetAccess=public,Hidden=true)
        DefaultLoss = @classreg.learning.loss.mincost;
        LabelPredictor = @classreg.learning.classif.ClassificationModel.minCost;
    end
    
    methods
        function cnames = get.ClassNames(this)
            cnames = labels(this.ClassSummary.ClassNames);
        end

        function cost = get.Cost(this)
            K = length(this.ClassSummary.ClassNames);
            cost = zeros(K);
            [~,pos] = ismember(this.ClassSummary.NonzeroProbClasses,...
                this.ClassSummary.ClassNames);
            cost(pos,pos) = this.ClassSummary.Cost;            
            unmatched = 1:K;
            unmatched(pos) = [];
            cost(:,unmatched) = NaN;
            cost(1:K+1:end) = 0;
        end

        function prior = get.Prior(this)
            K = length(this.ClassSummary.ClassNames);
            prior = zeros(1,K);
            [~,pos] = ismember(this.ClassSummary.NonzeroProbClasses,...
                this.ClassSummary.ClassNames);
            prior(pos) = this.ClassSummary.Prior;
        end
        
        function this = set.Prior(this,prior)
            this = setPrior(this,prior);
        end
        
        function this = set.Cost(this,cost)
            this = setCost(this,cost);
        end
        
        function st = get.ScoreTransform(this)
            st = classreg.learning.internal.convertScoreTransform(...
                this.PrivScoreTransform,'string',[]);
        end
        
        function this = set.ScoreTransform(this,st)
            this.PrivScoreTransform = ...
                classreg.learning.internal.convertScoreTransform(st,...
                'handle',numel(this.ClassSummary.ClassNames));
        end
    end
    
    methods(Access=protected,Abstract=true)
        s = score(this,X,varargin)
    end
    
    methods(Access=protected)
        function this = setPrior(this,~)
            error(message('stats:classreg:learning:classif:ClassificationModel:setPrior:Noop'));
        end
        
        function this = setCost(this,~)
            error(message('stats:classreg:learning:classif:ClassificationModel:setCost:Noop'));
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.Predictor(this,s);
            cnames = this.ClassNames;
            if ischar(cnames)
                s.ClassNames = cnames;
            else
                s.ClassNames = cnames';
            end
            s.ScoreTransform = this.ScoreTransform;
        end
                
        function [labels,posterior,cost] = predictEmptyX(this,X)
            D = numel(this.PredictorNames);
            if size(X,2)~=D
                error(message('stats:classreg:learning:classif:ClassificationModel:predictEmptyX:XSizeMismatch', D));
            end
            labels = repmat(this.ClassNames(1,:),0,1);
            K = numel(this.ClassSummary.ClassNames);
            posterior = NaN(0,K);
            cost = NaN(0,K);
        end        
    end

    methods(Hidden)
        function this = ClassificationModel(dataSummary,classSummary,scoreTransform)
            this = this@classreg.learning.Predictor(dataSummary);
            this.ClassSummary = classSummary;
            this.PrivScoreTransform = scoreTransform;
        end
        
        function this = setPrivatePrior(this,prior)
            if isempty(prior) || strncmpi(prior,'empirical',length(prior))
                error(message('stats:classreg:learning:classif:ClassificationModel:setPrivatePrior:EmpiricalOrEmptyPrior'));
            end
            this.ClassSummary.Prior = ...
                classreg.learning.classif.FullClassificationModel.processPrior(...
                prior,[],this.ClassSummary.ClassNames,...
                this.ClassSummary.NonzeroProbClasses);
        end
        
        function this = setPrivateCost(this,cost)
            this.ClassSummary.Cost = ...
                classreg.learning.classif.FullClassificationModel.processCost(...
                cost,this.ClassSummary.Prior,this.ClassSummary.ClassNames,...
                this.ClassSummary.NonzeroProbClasses);
        end
        
        function [X,C,W,Y,rowData] = prepareDataForLoss(this,X,Y,W,rowData,cleanRows)
            % Check X
            if ~isnumeric(X) || ~ismatrix(X)
                error(message('stats:classreg:learning:classif:ClassificationModel:prepareDataForLoss:BadXType'));
            end
            
            % Convert to class labels
            Y = classreg.learning.internal.ClassLabel(Y);
            
            % Check size
            N = size(X,1);
            if numel(Y)~=N
                error(message('stats:classreg:learning:classif:ClassificationModel:prepareDataForLoss:SizeXYMismatch'));
            end
            
            % Check weights
            if ~isfloat(W) || ~isvector(W) || length(W)~=N || any(W<0)
                error(message('stats:classreg:learning:classif:ClassificationModel:prepareDataForLoss:BadWeights', N));
            end
            W = W(:);
            
            % Clean rows for missing labels, zero prior and cost?
            if nargin<6 && N>0
                cleanRows = true;
            end
            
            % Any rowData present?
            if nargin>4 && ~isempty(rowData)
                haveRowData = true;
                if size(rowData,1)~=N
                    error(message('stats:classreg:learning:classif:ClassificationModel:prepareDataForLoss:SizeRowDataMismatch',N));
                end
            else
                haveRowData = false;
            end
            
            % Check for missing class labels
            t = ismissing(Y);
            if any(t) && cleanRows
                Y(t) = [];
                X(t,:) = [];
                W(t,:) = [];
                if haveRowData
                    rowData(t,:) = [];
                end
            end
            
            % Get a matrix of class counts
            C = classreg.learning.internal.classCount(this.ClassSummary.ClassNames,Y);
            
            % Remove observations for classes with zero probability
            zeroprior = this.Prior==0;
            if any(zeroprior) && cleanRows
                t = any(C(:,zeroprior),2);
                Y(t) = [];
                X(t,:) = [];
                C(t,:) = [];
                W(t,:) = [];
                if haveRowData
                    rowData(t,:) = [];
                end
            end
            
            % Remove observations for classes with zero cost
            zerocost = all(this.Cost==0,2)';
            if any(zerocost) && cleanRows
                t = any(C(:,zerocost),2);
                Y(t) = [];
                X(t,:) = [];
                C(t,:) = [];
                W(t,:) = [];
                if haveRowData
                    rowData(t,:) = [];
                end
            end
            
            % Normalize weights to the class prior
            if ~isempty(C)
                WC = bsxfun(@times,C,W);
                Wj = sum(WC,1);
                adjWFactor = zeros(1,numel(Wj));
                zeroprior = Wj==0;
                adjWFactor(~zeroprior) = this.Prior(~zeroprior)./Wj(~zeroprior);
                W = sum(bsxfun(@times,WC,adjWFactor),2);
            end
        end
    end
       
    methods
        function [labels,scores,cost] = predict(this,X,varargin)
        %PREDICT Predict response of the model.
        %   [LABEL,POSTERIOR,COST]=PREDICT(OBJ,X) returns predicted class labels
        %   LABEL, posterior probabilities POSTERIOR and misclassification costs
        %   COST for model OBJ and a matrix of predictors X. X must be a numeric
        %   matrix of size N-by-P, where P is the number of predictors used for
        %   training this model. Classification labels LABEL have the same type as
        %   Y used for training. Posterior probabilities POSTERIOR are an N-by-K
        %   numeric matrix for N observations and K classes. COST is an N-by-K
        %   matrix with predicted misclassification costs per class. The predicted
        %   label is assigned to the class with the minimal misclassification cost.
        %
        %   See also classreg.learning.classif.ClassificationModel.
  
            % Empty data
            if isempty(X)
                [labels,scores,cost] = predictEmptyX(this,X);
                return;
            end
            
            % Get scores from the compact class
            scores = score(this,X,varargin{:});
            
            % Transform scores and find the most probable class
            [labels,scores,cost] = this.LabelPredictor(this.ClassNames,...
                this.Prior,this.Cost,scores,this.PrivScoreTransform);
        end
    end
    
    methods
        function m = margin(this,X,Y,varargin)
        %MARGIN Classification margins.
        %   M=MARGIN(MODEL,X,Y) returns classification margins obtained by MODEL
        %   for matrix of predictors X and class labels Y. X must be a numeric
        %   matrix of size N-by-P, where P is the number of predictors used for
        %   training this model. Y must be of the same type as MODEL.Y and have N
        %   elements. Classification margin is the difference between
        %   classification score for the true class and maximal classification
        %   score for the false classes. The returned M is a numeric column-vector
        %   of length size(X,1).
        %
        %   See also classreg.learning.classif.ClassificationModel,
        %   classreg.learning.classif.ClassificationModel/predict.
        
            [X,C] = prepareDataForLoss(this,X,Y,ones(size(X,1),1),[],false);
            [~,Sfit] = predict(this,X,varargin{:});
            m = classreg.learning.loss.classifmargin(C,Sfit);
        end
        
        function e = edge(this,X,Y,varargin)
        %EDGE Classification edge.
        %   E=EDGE(MODEL,X,Y) returns classification edge obtained by MODEL for
        %   matrix of predictors X and class labels Y. X must be a numeric matrix
        %   of size N-by-P, where P is the number of predictors used for training
        %   this model. Y must be of the same type as MODEL.Y and have N elements.
        %   Classification edge is classification margin averaged over the entire
        %   data.
        %
        %   E=EDGE(OBJ,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'weights'   - Observation weights, a numeric vector of length
        %                     size(X,1). By default, all observation weights are
        %                     set to 1. If you supply weights, EDGE computes
        %                     weighted classification edge.
        %
        %   See also classreg.learning.classif.ClassificationModel,
        %   classreg.learning.classif.ClassificationModel/margin.
        
            % Get observation weights
            N = size(X,1);
            args = {'weights'};
            defs = {ones(N,1)};
            [W,~,extraArgs] = internal.stats.parseArgs(args,defs,varargin{:});
            
            % Prepare data
            [X,C,W] = prepareDataForLoss(this,X,Y,W);
            
            % Get observation margins for hypothesis H.
            [~,Sfit] = predict(this,X,extraArgs{:});

            % Check all arguments
            classreg.learning.internal.classifCheck(C,Sfit,W,[]);
            
            % Get edge            
            e = classreg.learning.loss.classifedge(C,Sfit,W,[]);
        end
        
        function l = loss(this,X,Y,varargin)
        %LOSS Classification error.
        %   L=LOSS(MODEL,X,Y) returns classification cost for model MODEL computed
        %   using matrix of predictors X and true class labels Y. X must be a
        %   numeric matrix of size N-by-P, where P is the number of predictors used
        %   for training this model. Y must be of the same type as MODEL.Y and have
        %   N elements.
        %
        %   L=LOSS(MODEL,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'lossfun'          - Function handle for loss, or string
        %                            representing a built-in loss function.
        %                            Available loss functions for classification:
        %                            'binodeviance', 'classiferror', 'exponential',
        %                            and 'mincost'. If you pass a function handle
        %                            FUN, LOSS calls it as shown below:
        %                                  FUN(C,S,W,COST)
        %                            where C is an N-by-K logical matrix for N rows
        %                            in X and K classes in the ClassNames property,
        %                            S is an N-by-K numeric matrix, W is a numeric
        %                            vector with N elements, and COST is a K-by-K
        %                            numeric matrix. C has one true per row for the
        %                            true class. S is a matrix of predicted scores
        %                            for classes with one row per observation,
        %                            similar to SCORE output from PREDICT. W is a
        %                            vector of observation weights. COST is a
        %                            matrix of misclassification costs. Default:
        %                            'mincost'
        %       'weights'          - Vector of observation weights. By default the
        %                            weight of every observation is set to 1. The
        %                            length of this vector must be equal to the
        %                            number of rows in X.
        %
        %   See also classreg.learning.classif.ClassificationModel,
        %   classreg.learning.classif.ClassificationModel/predict.  
            
            % Get loss function and observation weights
            N = size(X,1);
            args = {       'lossfun' 'weights'};
            defs = {this.DefaultLoss ones(N,1)};
            [funloss,W,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Prepare data
            [X,C,W] = prepareDataForLoss(this,X,Y,W);

            % Check loss function
            funloss = classreg.learning.internal.lossCheck(funloss,'classification');

            % Get observation margins for hypothesis H.
            [~,Sfit] = predict(this,X,extraArgs{:});

            % Check all arguments
            classreg.learning.internal.classifCheck(C,Sfit,W,this.Cost);

            % Get loss
            l = funloss(C,Sfit,W,this.Cost);
        end
    end
    
    methods(Static=true,Hidden=true)
        function [labels,scores,cost,classnum] = ...
                maxScore(classnames,Prior,Cost,scores,scoreTransform)
            scores = scoreTransform(scores);
            N = size(scores,1);
            notNaN = ~all(isnan(scores),2);
            [~,cls] = max(Prior);
            labels = repmat(classnames(cls,:),N,1);
            [~,classnum] = max(scores(notNaN,:),[],2);
            labels(notNaN,:) = classnames(classnum,:);
            cost = Cost(:,classnum)';
        end
        
        function [labels,scores,cost,classnum] = ...
                minCost(classnames,Prior,Cost,posterior,scoreTransform)
            cost = posterior*Cost;
            N = size(posterior,1);
            notNaN = ~all(isnan(cost),2);
            [~,cls] = max(Prior);
            labels = repmat(classnames(cls,:),N,1);
            [~,classnum] = min(cost(notNaN,:),[],2);
            labels(notNaN,:) = classnames(classnum,:);
            scores = scoreTransform(posterior);
        end
    end

end
