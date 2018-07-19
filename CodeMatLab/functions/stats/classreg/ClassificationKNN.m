
classdef ClassificationKNN < ...
        classreg.learning.classif.FullClassificationModel
%ClassificationKNN K Nearest Neighbors classification 
%   ClassificationKNN is a K Nearest Neighbors (KNN) classification model.
%   It can predict response for new data. It also stores data used for
%   training and can compute resubstitution predictions.
%
%   An object of this class cannot be created by calling the constructor.
%   Use ClassificationKNN.fit to create a ClassificationKNN object by
%   fitting the KNN classification model to training data.
%
%   ClassificationKNN properties:
%       NObservations         - Number of observations.
%       X                     - Matrix of predictors.
%       Y                     - True class labels.
%       W                     - Weights of observations.
%       ModelParams           - KNN classification model parameters.
%       PredictorNames        - Names of predictors.
%       ResponseName          - Name of the response variable.
%       ClassNames            - Names of classes in Y.
%       Cost                  - Misclassification costs.
%       Prior                 - Prior class probabilities.
%       ScoreTransform        - Transformation applied to predicted
%                               classification scores.
%       CategoricalPredictors - Indices of categorical predictors.
%       NSMethod              - Method for K nearest neighbors search.
%       NumNeighbors          - Number of nearest neighbors. 
%       Distance              - Distance metric. 
%       DistParameter         - Additional distance parameters.
%       IncludeTies           - Flag to include the tie neighbors.
%       BreakTies             - Method of breaking ties if more than one
%                               class has the same number of nearest points
%                               among the K nearest neighbors.
%  
%   ClassificationKNN methods:
%       fit (static)          - Fit a KNN classification model.
%       template (static)     - KNN template for FITENSEMBLE.
%       crossval              - Cross-validate this KNN classification model.
%       edge                  - Classification edge.
%       loss                  - Classification loss.
%       margin                - Classification margins.
%       predict               - Predicted response.
%       resubEdge             - Resubstitution classification edge.
%       resubLoss             - Resubstitution classification loss.
%       resubMargin           - Resubstitution classification margins.
%       resubPredict          - Resubstitution predicted response.  
%
%   Example: Grow a KNN classification model for Fisher's iris data.
%       load fisheriris
%       d = ClassificationKNN.fit(meas,species,'PredictorNames',...
%                    {'SL' 'SW' 'PL' 'PW'},'NumNeighbors',5)
                 

%  Copyright 2011 The MathWorks, Inc.

    properties(GetAccess=protected,SetAccess=protected)
        %NeighborSearcher object        
        NS = [];
        %saves the weight provided by the customer or come from the default
        %value. This weight is not normalized to the class prior.
        %
        PrivW = [];
        
                
    end
    
    
    properties(GetAccess=public,SetAccess=public,Dependent=true)
        %NumNeighbors Number of nearest neighbors. 
        %   A positive integer. The number of nearest neighbors for
        %   classifying each point when predicting.
        %
        %   See also ClassificationKNN
        NumNeighbors;
        
        %Distance Distance metric.
        %   A string specifying the built-in distance metric or a function
        %   handle.
        %
        %   See also ClassificationKNN
        Distance;
        
        %DistParameter Additional distance parameter.
        %   Value of Distance property Value
        % 
        %   'minkowski'                  A positive scalar indicating
        %                                the exponent of the minkowski     
        %                                distance. 
        %   'mahalanobis'                A positive definite matrix        
        %                                representing the covariance
        %                                matrix used for computing the
        %                                mahalanobis distance.
        %   'seuclidean'                 A vector representing the scale
        %                                value to use in computing the
        %                                'seuclidean' distance.
        %   otherwise                    Empty. 
        %
        %   See also ClassificationKNN
        DistParameter;
              
        %IncludeTies Whether to include the tie neighbors
        %   A logical value specifying whether to include all neighbors whose
        %   distance equal the Kth smallest distance.
        %
        %   See also ClassificationKNN
        IncludeTies;
        
        %DistanceWeight Distance weighting function.
        %   A string or a function handle specifying the distance weighting
        %   function. 
        %
        %   See also ClassificationKNN
        DistanceWeight;
                
        %BreakTies Method of breaking ties
        %   A string 'smallest', 'random' or 'nearest', specifying method of
        %   breaking ties if more than one class has the same smallest
        %   misclassification cost.
        %
        %   See also ClassificationKNN
        BreakTies;
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        % NSMethod Method for K nearest neighbors search
        %    A string 'kdtree' or 'exhaustive', specifying nearest neighbors
        %    search method.
        %
        %   See also ClassificationKNN
        NSMethod;
    end
    
    methods
               
        function nsmethod = get.NSMethod(this)
            nsmethod = this.ModelParams.NSMethod;
        end
        
                
        function dp = get.DistParameter(this)
            if ~isempty(this.ModelParams.Exponent)
                dp = this.ModelParams.Exponent;
            elseif ~isempty(this.ModelParams.Cov)
                dp = this.ModelParams.Cov;
            elseif ~isempty(this.ModelParams.Scale)
                dp = this.ModelParams.Scale;
            else
                dp =[];
            end
        end
        
        function dist = get.Distance(this)
            dist = this.ModelParams.Distance;
        end
        
        function this = set.Distance(this,distMetric)
           
            this.NS.Distance = distMetric;
            this.ModelParams.Distance = this.NS.Distance;
            %need to change the distance parameters
            i = find(strncmpi(this.ModelParams.Distance, ...
                {'minkowski','mahalanobis','seuclidean'},3));
            if isempty(i)
                 this.ModelParams.Exponent = [];
                 this.ModelParams.Cov =[];
                 this.ModelParams.Scale =[];
            elseif i == 1 %'min'
                 this.ModelParams.Exponent = this.NS.DistParameter;
                 this.ModelParams.Cov =[];
                 this.ModelParams.Scale =[];
            elseif i==2  %mah
                this.ModelParams.Cov= this.NS.DistParameter;
                this.ModelParams.Exponent = [];
                this.ModelParams.Scale = [];
            else %secu
                this.ModelParams.Scale= this.NS.DistParameter;
                this.ModelParams.Cov= [];
                this.ModelParams.Exponent = [];
            end
        end
           
        function this = set.DistParameter(this, para)
            this.NS.DistParameter = para;
            i = find(strncmpi(this.ModelParams.Distance, {'minkowski','mahalanobis','seuclidean'},3));
            if isempty(i)
                 error(message('stats:ClassificationKNN:set:DistParameter:InvalidDistanceParam'));
            elseif i == 1 %'min'
                 this.ModelParams.Exponent = para;
                 this.ModelParams.Cov =[];
                 this.ModelParams.Scale =[];
            elseif i == 2  %mah
                this.ModelParams.Cov= para;
                this.ModelParams.Exponent = [];
                this.ModelParams.Scale = [];
            else %secu
                this.ModelParams.Scale= para;
                this.ModelParams.Cov= [];
                this.ModelParams.Exponent = [];
            end
         end

        function K = get.NumNeighbors(this)
            K = this.ModelParams.NumNeighbors;
        end
        
        function this = set.NumNeighbors(this,K)
            if ~isscalar(K) || ~isnumeric(K) ||  ...
                    K <1 || K~=round(K)
                error(message('stats:ClassificationKNN:set:NumNeighbors:BadK'));
            end
            nx= size(this.X,1);
            this.ModelParams.NumNeighbors = min(K,nx);
        end
        
        function inTies = get.IncludeTies(this)
            inTies = this.ModelParams.IncludeTies;
        end
        
        function this = set.IncludeTies(this,tf)
            if ~islogical(tf) || ~isscalar(tf)
                error(message('stats:ClassificationKNN:set:IncludeTies:BadIncludeTies'));
            end
            this.ModelParams.IncludeTies = tf;
        end
        
        function inTies = get.BreakTies(this)
            inTies = this.ModelParams.BreakTies;
        end
        
        function this = set.BreakTies(this,breakties)
            if ~ischar(breakties)
                error(message('stats:ClassificationKNN:set:BreakTies:BadBreakTies'));
                   
            else
                breaktieList = {'smallest' 'nearest','random'};
                i = find(strncmpi(breakties,breaktieList,length(breakties)));
                if isempty(i)
                    error(message('stats:ClassificationKNN:set:BreakTies:BadBreakTies'));
                else
                    this.ModelParams.BreakTies = breaktieList{i};
                end
            end
        end
        
        function distWgt = get.DistanceWeight(this)
            distWgt = this.ModelParams.DistanceWeight;
        end
        
        function this = set.DistanceWeight(this,wgt)
            if ischar(wgt)
                wgtList = {'equal','inverse','squaredinverse'};
                i = find(strncmpi(wgt, wgtList,length(wgt)));
                if isempty(i)
                    error(message('stats:ClassificationKNN:set:DistanceWeight:BadDistanceWeight'));
                else
                    this.ModelParams.DistanceWeight = wgtList{i};
                end
              
            elseif isa (wgt,  'function_handle')
                this.ModelParams.DistanceWeight = wgt;
            else
                error(message('stats:ClassificationKNN:set:DistanceWeight:BadDistanceWeight'));           
            end          
        end
      
    end
    
    methods(Access=protected)
           %This function normalzied observation weights to be consistent with the class priors 
       %i.e., the summation of observation weights in each class will equal
       %to the class prior
       function this = normalizeWeights(this)
            C = classreg.learning.internal.classCount(...
                this.ClassSummary.NonzeroProbClasses,this.PrivY);
         
            WC = bsxfun(@times,C,this.PrivW);
            Wj = sum(WC,1);
            this.W = sum(bsxfun(@times,WC,this.ClassSummary.Prior./Wj),2);
       end
         
       function this = setPrior(this,prior)
             % call setPrivatePrior in 
             % +classreg\+learning\+classif\ClassificationModel.m
            this = setPrivatePrior(this,prior);  
%             Each time when the prior is set, we need to normalized the
%             observation weights using the new prior
            this = normalizeWeights(this);
        end
         
         function this = setCost(this,cost)
             % call setPrivateCost in 
             % +classreg\+learning\+classif\ClassificationModel.m
             this = setPrivateCost(this,cost);       
        end

       %add properties for display
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.FullClassificationModel(this,s);
            s.Distance = this.Distance;
            s.NumNeighbors = this.NumNeighbors;  
        end
        
        %must return a N by K(#of classes) matrix for N observation and
        %K classess
        function [s,gindex,CIDX]=score(this,X,varargin)
                   
 
            W = this.W;
            NG = length(this.ClassSummary.ClassNames);
           
            includeTies = this.ModelParams.IncludeTies;
            gindex = grp2idx(this.PrivY,this.ClassSummary.ClassNames); 
            distWgtList ={'equal','inverse','squaredinverse'};
            distanceWeightFun = this.ModelParams.DistanceWeight;
            distWgtIdx = find(strncmpi(distanceWeightFun,distWgtList,3));
           
            [CIDX,dist] = knnsearch(this.NS, X,'k',this.ModelParams.NumNeighbors,...
                       'includeTies', includeTies);
            
            NX= size(X,1);
          
            % Returns scores for all classes but compute the scores only
            % for classes with observations in the training data
          
            %  count is a NX-by-NG matrix saving the count of the neighbors
            %  in each class for each query point, NY is the number of points
            %  in test and NG is number of groups
            count =zeros(NX,NG);
            if (includeTies)
                count = zeros(NG,NX);
                if isa(distanceWeightFun,'function_handle')
                    
                    try
                        distWgt = feval(this.ModelParams.DistanceWeight,...
                            dist{1});
                        
                    catch ME
                        if strcmp('MATLAB:UndefinedFunction', ME.identifier) ...
                                && ~isempty(strfind(ME.message, func2str(distanceWeightFun)))
                            error(message('stats:ClassificationKNN:Score:DistanceFunctionNotFound', func2str( distanceWeightFun )));
                        end
                        
                    end
                    if ~isnumeric(distWgt)
                        error(message('stats:ClassificationKNN:Score:OutputBadType'));
                    end
                   
                    for outer =1:NX
                        %nan values in the dist will be ignored
                        numNeighbors = sum(~isnan(dist{outer}));
                        tempCIDX =CIDX{outer}(1:numNeighbors);
                        tempIDX =gindex(tempCIDX);
                        tempDist = dist{outer}(1:numNeighbors);
                        obsWgt = W(tempCIDX);              
                        distWgt = feval(this.ModelParams.DistanceWeight,tempDist );
                        if (any(distWgt<0))
                             error(message('stats:ClassificationKNN:Score:NegativeDistanceWgt'));
                        end
                        wgt = obsWgt .* distWgt';
                        wgt(isnan(wgt)) = 0; %don't consider neighbors with NaN weights
                        %Both tempIDX and wgt are numneighbor-by-one vectors 
                        count(:,outer) = ...
                                accumarray(tempIDX,wgt,[NG,1]);              
                    end
         
                elseif distWgtIdx == 1 %equal weight
                    
                    for outer =1:NX
                        %nan values in the dist will be ignored
                        numNeighbors = sum(~isnan(dist{outer}));
                        tempCIDX =CIDX{outer}(1:numNeighbors);
                        tempIDX =gindex(tempCIDX);
                        wgt = W(tempCIDX);
                        count(:,outer) = ...
                                accumarray(tempIDX,wgt,[NG,1]);                         
                    end
                       
                else % inverse weight or squared inverse weight
                   
                    if distWgtIdx==2 %'inverse'
                        e = 1;
                    else
                        e = 2;
                    end
                    for outer =1:NX
                        %nan values in the dist will be ignored
                        numNeighbors = sum(~isnan(dist{outer}));
                        tempCIDX =CIDX{outer}(1:numNeighbors);
                        tempIDX =gindex(tempCIDX);
                        
                        obsWgt = W(tempCIDX);
                        tempDist = dist{outer}(1:numNeighbors);
                        distWgt = wgtFunc(tempDist,e);
                        wgt = obsWgt .* distWgt';
                        count(:,outer) = ...
                                accumarray(tempIDX,wgt,[NG,1]);

                    end  
                   
                end
                count = count';            
            else % includeTie is false
                                   
                %Find the number of valid neighbors
                %The columns in matrix dist with nan values in the whole column
                %will be ignored. 
                %note that some distance, (e.g., cosine) may have NaN
                %values even the data doesn't have NaNs. Therefor the
                %following step is needed even we remove the observations
                %with NaN values in the training data.
                numNeighbors = sum((~all(isnan(dist),1)));
                if numNeighbors > 0
                    dist(:,numNeighbors+1:end)=[]; %don't consider invalid neighbors
                    CIDX(:,numNeighbors+1:end)=[];
                    
                    if isa(distanceWeightFun,'function_handle')
                        try
                            distWgt = feval(this.ModelParams.DistanceWeight,dist(1,:));
                        catch ME
                            if strcmp('MATLAB:UndefinedFunction', ME.identifier) ...
                                    && ~isempty(strfind(ME.message, func2str(distanceWeightFun)))
                                error(message('stats:ClassificationKNN:Score:DistanceFunctionNotFound',...
                                    func2str( distanceWeightFun )));
                            end
                        end
                        if ~isnumeric(distWgt)
                            error(message('stats:ClassificationKNN:Score:OutputBadType'));
                        end
                        distWgt = feval(this.ModelParams.DistanceWeight,dist );
                        if any(distWgt(:)< 0)
                             error(message('stats:ClassificationKNN:Score:NegativeDistanceWgt'));
                       end
                    elseif distWgtIdx==1 %equal weight
                        distWgt = ones(NX,numNeighbors);
                        %weights for neighbors which has NaN distance to
                        %the test point is set to zero,so that this
                        %neighbr will be ingored.
                        distWgt(isnan(dist)) = 0;
                    elseif distWgtIdx==2 %inverse
                        distWgt = wgtFunc(dist,1);
                    else
                        distWgt = wgtFunc(dist,2);
                    end
                    
                    %CNeighbor is a matrix with size NX by numNeighbors
                    %representing the class labels for nearest neighbors             
                    CNeighbor = gindex(CIDX);
                    obsWgt= W(CIDX);
                    
                    if (NX==1) && numNeighbors > 1
                        CNeighbor = CNeighbor';
                        obsWgt = obsWgt';
                    end
                    
                    %wgt is a matrix with size NX by numNeighbors
                    wgt = distWgt .* obsWgt;
                    %don't consider neighbors with NaN weights
                    %This NaN values may happen, for example, if the number
                    %of valid training points is less than the number of
                    %requsted neighbors.
                    wgt(isnan(wgt)) = 0; 
                    if numNeighbors > 5
                        %when the number of neighbors is bigger than 5,
                        %using accumarray is faster than using the loop
                        count =zeros(NG,NX);
                        wgt=wgt';
                        CNeighbor = CNeighbor';
                        for i = 1:NX
                            count(:,i) = ...
                                accumarray(CNeighbor(:,i),wgt(:,i),[NG,1]);
                        end
                        count = count';
                    else %
                        count = zeros(NX,NG);
                        for outer = 1:NX
                            for inner = 1:numNeighbors
                                count(outer,CNeighbor(outer,inner))=...
                                    count(outer,CNeighbor(outer,inner))+wgt(outer,inner);
                            end
                        end
                    end
                end
            end
         
            if isa(distanceWeightFun,'function_handle')
                %deal with inf weight values when distance weight
                %function is function handle
                infCountRow = any(isinf(count),2);
                s = count;
                s(infCountRow,:) = 0;
                s(isinf(count)) = 1;
                %to deal with the case that one point has two
                %points with infinite weight but coming from
                %different classes.
                s=bsxfun(@rdivide, s, sum(s,2));
            else
                s=bsxfun(@rdivide, count, sum(count,2));
            end
         
        end      
        
    end
    
    methods
        function partModel = crossval(this,varargin)
        %CROSSVAL Cross-validate this model.
        %   CVMODEL=CROSSVAL(MODEL) builds a partitioned model CVMODEL from model
        %   MODEL represented by a full object for classification. You can then
        %   assess the predictive performance of this model on cross-validated data
        %   using methods and properties of CVMODEL. By default, CVMODEL is built
        %   using 10-fold cross-validation on the training data. CVMODEL is of
        %   class ClassificationPartitionedModel.
        %
        %   CVMODEL=CROSSVAL(MODEL,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %      'kfold'      - Number of folds for cross-validation, a numeric
        %                     positive scalar; 10 by default.
        %      'holdout'    - Holdout validation uses the specified
        %                     fraction of the data for test, and uses the rest of
        %                     the data for training. Specify a numeric scalar
        %                     between 0 and 1.
        %      'leaveout'   - If 'on', use leave-one-out cross-validation.
        %      'cvpartition' - An object of class CVPARTITION; empty by default. If
        %                      a CVPARTITION object is supplied, it is used for
        %                      splitting the data into subsets.
        %
        %   See also classreg.learning.classif.FullClassificationModel,
        %   cvpartition,
        %   classreg.learning.partition.ClassificationPartitionedModel.

     
            idxBaseArg = find(ismember(lower(varargin(1:2:end)),...
                classreg.learning.FitTemplate.AllowedBaseFitObjectArgs));
            if ~isempty(idxBaseArg)
                error(message('stats:classreg:learning:classif:FullClassificationModel:crossval:NoBaseArgs', varargin{ 2*idxBaseArg - 1 }));
            end
            temp = classreg.learning.FitTemplate.make(this.ModelParams.Method,...
                'type','classification','scoretransform',this.PrivScoreTransform,...
                'modelparams',this.ModelParams,'crossval','on',varargin{:});
          %  ClassificationKNN/crossval must pass prior but not cost o the
          %  constructor of the partitioned model. It must assign cost
          %  after the partitioned model is constructed. 
            partModel = fit(temp,this.X,this.Y,'weights',this.W,...
                'predictornames',this.DataSummary.PredictorNames,...
                'categoricalpredictors',this.CategoricalPredictors,...
                'responsename',this.ResponseName,...
                'classnames',this.ClassNames,'prior',this.Prior);
            partModel.Cost = this.Cost;
        end
        
        function [label,posteriors,cost] = predict(this,X)
        %PREDICT Predict response of KNN classification model.
            %   [LABEL,POSTERIOR,COST]=PREDICT(KNN,X) returns predicted
            %   class labels LABEL, posterior probabilities POSTERIOR and
            %   misclassification costs COST for KNN and a matrix of
            %   predictors X. X must be a numeric matrix of size N-by-P,
            %   where P is the number of predictors used for training this
            %   model. Classification labels LABEL have the same type as Y
            %   used for training. Posterior probabilities POSTERIOR are an
            %   N-by-K numeric matrix for N observations and K classes.
            %   COST is an N-by-K matrix with predicted misclassification
            %   costs per class. The predicted label is assigned to the
            %   class with the minimal misclassification cost.
            %
            %   See also ClassificationKNN
            
            %override the predict method in ClassificationModel to add the
            %third output COST
            % Empty data
            if isempty(X)
                [label,posteriors] = predictEmptyX(this,X);
                cost = NaN(0,numel(this.ClassSummary.ClassNames));
                return;
            end
            
            % Get posterior probabilities(P(x|k)) from score() method
             breakTieFlag= find(strncmpi(this.BreakTies,{'random','nearest'},3));
              %if 'BreakTies' is 'smallest' or 'random' or only 1 nearest
              %neighbor is used
              if isempty(breakTieFlag) || breakTieFlag == 1 ||...
                      this.NumNeighbors ==  1
                 posteriors = score(this,X);
              else
                 [posteriors,gindex,CIDX] = score(this,X);
              end
             % Transform posterior, compute expected cost and find the class
             % with minimal predicted cost
%             [label,posterior,cost] = this.LabelPredictor(this.ClassNames,...
%                 this.Prior,this.Cost,posterior,this.PrivScoreTransform);
            cost = posteriors*this.Cost;
            N = size(posteriors,1);
            notNaN = ~all(isnan(cost),2);
            [~,cls] = max(this.Prior);
            label = repmat(this.ClassNames(cls,:),N,1);
            minCost = nan(N,1);
            [minCost(notNaN),classNum] = min(cost(notNaN,:),[],2);
            label(notNaN,:) = this.ClassNames(classNum,:);
            posteriors = this.PrivScoreTransform(posteriors);
            
            %deal with BreakTies
            if ~isempty(breakTieFlag) && this.NumNeighbors > 1
                
                notNanRows = find(notNaN);
                if breakTieFlag==1 % BreakTies is 'random'
                    for i =1:numel(notNanRows)
                        ties = cost(notNanRows(i),:) == minCost(notNanRows(i));
                        numTies = sum(ties);
                        if  numTies > 1 % there existed ties
                            choice = find(ties);
                            tb = randsample(numTies,1);
                            label(notNanRows(i)) = choice(tb);
                        end
                    end
                else % BreakTies is 'nearest'
                    
                    if ~this.IncludeTies  %this.IncludeTies is false
                        CNeighbor = gindex(CIDX);
                        for i =1:numel(notNanRows)
                            ties = cost(notNanRows(i),:) == minCost(notNanRows(i));
                            numTies = sum(ties);
                            if  numTies > 1 % there existed ties
                                choice = find(ties);
                                for inner = 1:this.NumNeighbors
                                    if ismember(CNeighbor(notNanRows(i),inner),choice)
                                        label(notNanRows(i)) = CNeighbor(notNanRows(i),inner);
                                        break
                                        
                                    end
                                end
                            end
                        end %
                    else %this.IncludeTies is true
                        
                        for i =1:numel(notNanRows)
                            ties = cost(notNanRows(i),:) == minCost(notNanRows(i));
                            numTies = sum(ties);
                            if  numTies > 1 % there existed ties
                                choice = find(ties);
                                for inner = 1:this.NumNeighbors
                                    tempCNeighbor =gindex(CIDX{notNanRows(i)});
                                    if ismember(tempCNeighbor(inner),choice)
                                        label(notNanRows(i)) = tempCNeighbor(inner);
                                        break
                                        
                                    end
                                end
                            end
                        end %if i==1
                    end
                end
            end
        end %~isempty(i)
                           
    end
    
    methods(Hidden)
        %constructor
        function this = ClassificationKNN(X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform)
            % Protect against being called by the user
            if nargin~=7 || ischar(W)
                error(message('stats:ClassificationKNN:ClassificationKNN:DoNotUseConstructor'));
            end
            
            [nx,nDims] = size(X);
            if (modelParams.NumNeighbors > nx)
                modelParams.NumNeighbors = nx;
            end          
            if ~ (isempty(dataSummary.CategoricalPredictors) || ...
                    (length(dataSummary.CategoricalPredictors) == nDims ...
                    && all(dataSummary.CategoricalPredictors==(1: nDims))))
                error(message('stats:ClassificationKNN:ClassificationKNN:BadCategoricalPre')); 
            end
            
            % Base constructor
            this =this@classreg.learning.classif.FullClassificationModel(...
                X,Y,W,modelParams,dataSummary,classSummary,scoreTransform);
            %we don't have class CompactClassificationDiscriminant
            %             this = this@classreg.learning.classif.CompactClassificationDiscriminant(...
            %                 dataSummary,classSummary,scoreTransform,[]);
            %this.LabelPredictor = @classreg.learning.classif.ClassificationModel.minCost;
            %Train the model on this.X, this.Y
           
            nsmethod= this.ModelParams.NSMethod;
            distance = this.ModelParams.Distance;
            p = this.ModelParams.Exponent;
            cov = this.ModelParams.Cov;
            scale = this.ModelParams.Scale;
            bucket = this.ModelParams.BucketSize;
            if strncmpi(nsmethod,'kdtree',length(nsmethod))
                this.NS = KDTreeSearcher(this.X, ...
                    'distance',distance, 'p',p,'BucketSize',bucket);
                this.ModelParams.NSMethod = 'kdtree';
                this.ModelParams.Scale = [];
                this.ModelParams.Cov = [];
                this.ModelParams.Exponent = this.NS.DistParameter;
            else
                %checkNegativeDistance is set to true to disallow negative distance
                %values.
                this.NS = ExhaustiveSearcher(this.X,...
                    'distance',distance, 'p',p,'cov',cov,...
                    'scale',scale,'checkNegativeDistance',true);
                this.ModelParams.NSMethod = 'exhaustive';
                this.ModelParams.BucketSize = [];
            end
            
            this.ModelParams.Distance = this.NS.Distance;
            
            %saves the observation weights that comes from the input,
            %This PrivW is not normalized to be consistent with class priors.
            %This is to prevent losing some observation weights when some class
            %priors are setting to zero.
            this.PrivW = this.W; 

        end
        
    end
    
    methods(Static)
        function this = fit(X,Y,varargin)
% FIT fit KNN classification model
%   KNN = CLASSIFICATIONKNN.FIT(X,Y) returns a KNN classification  model
%   for predictors X and class labels Y.
%
% X must be an N-by-P matrix of predictors with one row per observation and
% one column per predictor. Y must be an array of N class labels. Y can be
% a categorical array (nominal or ordinal), character array, logical
% vector, numeric vector, or cell array of strings. If Y is a character
% array, it must have one class label per row.
%
% KNN is a KNN classification model. If you use one of the following five
% options, KNN is of class ClassificationPartitionedModel: 'crossval',
% 'kfold', 'holdout', 'leaveout' or 'cvpartition'. Otherwise, KNN is of
% class ClassificationKNN.
%
%   KNN=CLASSIFICATIONKNN.FIT(X,Y,'PARAM1',val1,'PARAM2',val2,...)
%   specifies optional parameter name/value pairs:
%       'CategoricalPredictors' - List of categorical predictors. Pass
%                                 'CategoricalPredictors' as [] or 'all'.
%                                 Use [] to indicate none of predictors are
%                                 categorical. Use 'all' to indicate all
%                                 the predictors are categorical.
%                                 Default: []
%       'ClassNames'            - Array of class names. Use the
%                                 data type that exists in Y. You can use
%                                 this argument to order the classes or
%                                 select a subset of classes for training.
%                                 Default: The class names that exist in Y.
%       'Cost'                  - Square matrix, where COST(I,J) is the
%                                 cost of classifying a point into class J
%                                 if its true class is I.
%                                 Alternatively, COST can be a structure S
%                                 having two fields: S.ClassificationCosts
%                                 containing the cost matrix C, and
%                                 S.ClassNames containing the class names
%                                 and defining the ordering of classes used
%                                 for the rows and columns of the cost
%                                 matrix. For S.ClassNames use the data
%                                 type that exists in Y. As an alternative,
%                                 you can assign to the Cost property of
%                                 KNN. Default: COST(I,J)=1 if I~=J, and
%                                 COST(I,J)=0 if I=J.
%       'CrossVal'              - If 'on', creates a cross-validated KNN
%                                 with 10 folds. You can use 'kfold',
%                                 'holdout', 'leaveout' and 'cvpartition'
%                                 parameters to override this
%                                 cross-validation setting. You can only
%                                 use one of these four options ('kfold',
%                                 'holdout', 'leaveout' and 'cvpartition')
%                                 at a time when creating a cross-validated
%                                 model. As an alternative, you can
%                                 cross-validate later using CROSSVAL
%                                 method for ClassificationKNN. Default:
%                                 'off'
%       'CVPartition'           - A partition created with  CVPARTITION to
%                                 use in cross-validated KNN model.
%       'Holdout'               - Holdout validation uses the specified
%                                 fraction of the data for test, and uses
%                                 the rest of the data for training.
%                                 Specify a numeric scalar between 0 and 1.
%       'KFold'                 - Number of folds to use in cross-validated
%                                 KNN classification model, a positive
%                                 integer. Default: 10
%       'Leaveout'              - Use leave-one-out cross-validation by
%                                 setting to 'on'.
%       'Prior'                 - Prior probabilities for each class.
%                                 Specify as one of:
%                                   * A string:
%                                     - 'empirical' determines class
%                                       probabilities from class
%                                       frequencies in Y
%                                     - 'uniform' sets all class
%                                       probabilities equal
%                                   * A vector (one scalar value for each
%                                     class)
%                                   * A structure S with two fields:
%                                     S.ClassProbs containing a vector of
%                                     class probabilities, and S.ClassNames
%                                     containing the class names and
%                                     defining the ordering of classes used
%                                     for the elements of this vector.
%                                 If you pass numeric values,
%                                 ClassificationKNN.fit normalizes
%                                 them to add up to one. As an
%                                 alternative, you can assign to
%                                 the Prior property. Default:
%                                 'empirical'
%       'ResponseName'          - Name of the response variable Y, a
%                                 string. Default: 'Y'
%       'ScoreTransform'        - Function handle for transforming scores,
%                                 or string representing a built-in
%                                 transformation function. Available
%                                 functions: 'symmetric', 'invlogit',
%                                 'ismax', 'symmetricismax', 'none',
%                                 'logit', 'doublelogit', 'symmetriclogit',
%                                 and 'sign'.
%                                 Default: 'none'
%       'Weights'               - Vector of observation weights, one weight
%                                 per observation. Default:
%                                 ones(size(X,1),1)
%       'NumNeighbors'          - A positive integer specifying the number
%                                 of nearest neighbors in X for classifying
%                                 each point when predicting. Default: 1.
%       'NSMethod'              - Nearest neighbors search method.
%                                 Value is either:
%                                   - 'kdtree' uses a kd-tree to find
%                                     nearest neighbors. 'kdtree' is only
%                                     valid when the distance metric
%                                     is one of the following metrics:
%                                       - 'euclidean'
%                                       - 'cityblock'
%                                       - 'minkowski'
%                                       - 'chebyshev'
%                                   - 'exhaustive' uses the exhaustive
%                                     search algorithm. The distance values
%                                     from all the points in X to each
%                                     point in Y are computed to find
%                                     nearest neighbors.
%                                   Default is 'kdtree' when the number of
%                                   columns of X is not greater  than 10, X
%                                   is not sparse, and the distance metric
%                                   is one of the above 4 metrics;
%                                   otherwise, default is 'exhaustive'.
%       'IncludeTies'           - A logical value. Use true to include all
%                                 neighbors whose distance equal the Kth
%                                 smallest distance. Use false to include
%                                 exactly K nearest neighbors.
%       'DistanceWeight'        - A string or a function handle specifying
%                                 the distance weighting function. The
%                                 choices of the string are:
%                                   - 'equal': Each neighbor gets equal
%                                      weight (default).
%                                   - 'inverse': Each neighbor gets weight
%                                      1/d, where d is the distance between
%                                      this neighbor and the point being
%                                      classified.
%                                   - 'squaredinverse': Each neighbor gets
%                                     weight 1/d^2, where d is the distance
%                                     between this neighbor and the point
%                                     being classified.
%                                   - A distance weighting function
%                                     specified using @. A distance
%                                     weighting function must be of the
%                                     form:
%
%                                     function DW = DISTWGT(D)
%
%                                     taking as argument a matrix D and
%                                     returning a matrix of distance weight
%                                     DW. D and DW can only contains
%                                     non-negative numerical values. DW must
%                                     have the same size as D. DW(I,J) is
%                                     the weight computed based on D(I,J).
%        'BreakTies'            - Method of breaking ties if more than one
%                                 class has the same smallest cost. Choices
%                                 are:
%                                   - 'smallest': Assign the point to the
%                                     class with the smallest index. This
%                                     is default.
%                                   - 'nearest': Assign the point to the
%                                     class of its nearest neighbor.
%                                   - 'random': Randomly pick a class
%                                      out the classes with the smallest
%                                      cost.
%        'Distance'              - A string or a function handle specifying
%                                  the distance metric. The value can be
%                                  one of the following:
%                                  'euclidean'   
%                                     - Euclidean distance. This is default
%                                       if 'CategoricalPredictors' is
%                                       empty.
%                                  'seuclidean'
%                                     - Standardized Euclidean distance.
%                                       Each coordinate difference between
%                                       X and a query point is divided by
%                                       an element of  vector S. The
%                                       default value of is the standard
%                                       deviation computed from X,
%                                       S=NANSTD(X). To specify another
%                                       value for S, use the 'Scale'
%                                       argument.
%                                  'cityblock'  
%                                     - City Block distance.
%                                  'chebychev'  
%                                     - Chebychev distance (maximum
%                                       coordinate  difference).
%                                  'minkowski'  
%                                     - Minkowski distance. The default
%                                       exponent is 2. To specify a
%                                       different exponent, use the 'P'
%                                       argument.
%                                  'mahalanobis' 
%                                     - Mahalanobis distance, computed
%                                       using a positive definite
%                                       covariance matrix C. The default
%                                       value of C is the sample covariance
%                                       matrix of X, as computed by
%                                       NANCOV(X). To specify another value
%                                       for C, use the 'Cov' argument.
%                                  'cosine' 
%                                     - One minus the cosine of the
%                                       included angle between observations
%                                       (treated as vectors).
%                                   'correlation' 
%                                     - One minus the sample linear
%                                       correlation between observations
%                                       (treated as sequences of values).
%                                   'spearman'   
%                                     - One minus the sample Spearman's
%                                       rank correlation between
%                                       observations (treated as sequences
%                                       of values).
%                                   'hamming'     
%                                     - Hamming distance, percentage of
%                                       coordinates that differ. This is
%                                       default if 'CategoricalPredictors'
%                                       is set to 'all'.
%                                   'jaccard'     
%                                     - One minus the Jaccard coefficient,
%                                       the percentage of nonzero
%                                       coordinates that differ.
%                                   function     
%                                     - A distance function specified using
%                                       @ (for example @DISTFUN). A
%                                       distance function must be of the
%                                       form
%
%                                       function D2 = DISTFUN(ZI, ZJ),
%
%                                       taking as arguments a 1-by-N vector
%                                       ZI containing a single row of X or
%                                       Y, an M2-by-N matrix ZJ containing
%                                       multiple rows of X or Y, and
%                                       returning an M2-by-1 vector of
%                                       distances D2, whose Jth element is
%                                       the distance between the
%                                       observations ZI and ZJ(J,:).
%       'Exponent'             - A positive scalar indicating the
%                                   exponent of Minkowski distance. This
%                                   argument is only valid when
%                                  'Distance' is 'minkowski'. Default: 2.
%       'Cov'                   - A positive definite matrix indicating the
%                                 covariance matrix when computing the
%                                 Mahalanobis distance. This argument is
%                                 only valid when 'Distance' is
%                                 'mahalanobis'. Default is NANCOV(X).
%       'Scale'                 - A vector S containing non-negative
%                                 values, with length equal to the number
%                                 of columns in X. Each  coordinate
%                                 difference between X and a query point is
%                                 divided by the corresponding element of
%                                 S. This argument is only valid when
%                                 'Distance' is 'seuclidean'. Default is
%                                 NANSTD(X).
%       'BucketSize'            - The maximum number of data points in the
%                                 leaf node of the kd-tree (default is 50).
%                                 This argument is only meaningful when
%                                 'NSMethod' is 'kdtree'.
%
%    See also ClassificationKNN

            Nfold = classreg.learning.generator.Partitioner.processArgs(varargin{:});
            docv = ~isempty(Nfold);
            
            % Set prior to an empty array. For a single KNN classifier, it
            % can be filled with something.
            prior = [];
            
            if ~docv %single KNN
                % For a single KNN, exclude cost and prior from input arguments.
                % If we pass a zero prior for a class, PrivW for this class
                % will be set to 0 in the ClassificationKNN constructor.
                % To prevent this, we pass a default empirical prior to fit
                % and reset the prior after fitting.
                args = {'prior' 'cost'};
                defs = {     []  []};
                [prior,cost,~,fitArgs] = internal.stats.parseArgs(args,defs,varargin{:});
              
            else % cross-validated KNN
                % For cross-validated KNN, we pass the prior specified by
                % the customer because a cross-validated KNN does not allow
                % assignment to Prior property.
                 args = {'cost'};
                 defs = {    []};
                [cost,~,fitArgs] = internal.stats.parseArgs(args,defs,varargin{:});
            end
            
            % Fit
            temp = classreg.learning.FitTemplate.make(...
                'KNN','type','classification',fitArgs{:});
            this = fit(temp,X,Y);
            
            % Assign prior and cost
            if ~isempty(prior) && ~strncmpi(prior,'empirical',length(prior))
                %this.W will be normalized based on class priors 
                this.Prior = prior;
            end
            if ~isempty(cost)
                this.Cost = cost;
           end
        end
        
         function temp = template(varargin)
%TEMPLATE Create a classification template.
%   T=CLASSIFICATIONKNN.TEMPLATE() returns a KNN classification
%   template suitable for use in the FITENSEMBLE function.
%
%   T=CLASSIFICATIONKNN.TEMPLATE('PARAM1',val1,'PARAM2',val2,...)
%   specifies optional parameter name/value pairs:
%       'NumNeighbors'          - A positive integer specifying the number
%                                 of nearest neighbors in X for classifying
%                                 each point when predicting. Default: 1.
%       'NSMethod'              - Nearest neighbors search method.
%                                 Value is either:
%                                   - 'kdtree' uses a kd-tree to find
%                                     nearest neighbors. 'kdtree' is only
%                                     valid when the distance metric
%                                     is one of the following metrics:
%                                       - 'euclidean'
%                                       - 'cityblock'
%                                       - 'minkowski'
%                                       - 'chebyshev'
%                                   - 'exhaustive' uses the exhaustive
%                                     search algorithm. The distance values
%                                     from all the points in X to each
%                                     point in Y are computed to find
%                                     nearest neighbors.
%                                   Default is 'kdtree' when the number of
%                                   columns of X is not greater  than 10, X
%                                   is not sparse, and the distance metric
%                                   is one of the above 4 metrics;
%                                   otherwise, default is 'exhaustive'.
%       'IncludeTies'           - A logical value. Use true to include all
%                                 neighbors whose distance equal the Kth
%                                 smallest distance. Use false to include
%                                 exactly K nearest neighbors.
%       'DistanceWeight'        - A string or a function handle specifying
%                                 the distance weighting function. The
%                                 choices of the string are:
%                                   - 'equal': Each neighbor gets equal
%                                      weight (default).
%                                   - 'inverse': Each neighbor gets weight
%                                      1/d, where d is the distance between
%                                      this neighbor and the point being
%                                      classified.
%                                   - 'squaredinverse': Each neighbor gets
%                                     weight 1/d^2, where d is the distance
%                                     between this neighbor and the point
%                                     being classified.
%                                   - A distance weighting function
%                                     specified using @. A distance
%                                     weighting function must be of the
%                                     form:
%
%                                     function DW = DISTWGT(D)
%
%                                     taking as argument a matrix D and
%                                     returning a matrix of distance weight
%                                     DW. D and DW can only contain
%                                     non-negative numerical values. DW must
%                                     have the same size as D. DW(I,J) is
%                                     the weight computed based on D(I,J).
%        'BreakTies'            - Method of breaking ties if more than one
%                                 class has the same smallest cost. Choices
%                                 are:
%                                   - 'smallest': Assign the point to the
%                                     class with the smallest index. This
%                                     is default.
%                                   - 'nearest': Assign the point to the
%                                     class of its nearest neighbor.
%                                   -  'random': Randomly pick a class
%                                      out the classes with the smallest
%                                      cost.
%        'Distance'              - A string or a function handle specifying
%                                  the distance metric. The value can be
%                                  one of the following:
%                                  'euclidean'   
%                                     - Euclidean distance (default).
%                                  'seuclidean'
%                                     - Standardized Euclidean distance.
%                                       Each coordinate difference between
%                                       X and a query point is divided by
%                                       an element of  vector S. The
%                                       default value of is the standard
%                                       deviation computed from X,
%                                       S=NANSTD(X). To specify another
%                                       value for S, use the 'Scale'
%                                       argument.
%                                  'cityblock'  
%                                     - City Block distance.
%                                  'chebychev'  
%                                     - Chebychev distance (maximum
%                                       coordinate  difference).
%                                  'minkowski'  
%                                     - Minkowski distance. The default
%                                       exponent is 2. To specify a
%                                       different exponent, use the 'P'
%                                       argument.
%                                  'mahalanobis' 
%                                     - Mahalanobis distance, computed
%                                       using a positive definite
%                                       covariance matrix C. The default
%                                       value of C is the sample covariance
%                                       matrix of X, as computed by
%                                       NANCOV(X). To specify another value
%                                       for C, use the 'Cov' argument.
%                                  'cosine' 
%                                     - One minus the cosine of the
%                                       included angle between observations
%                                       (treated as vectors).
%                                   'correlation' 
%                                     - One minus the sample linear
%                                       correlation between observations
%                                       (treated as sequences of values).
%                                   'spearman'   
%                                     - One minus the sample Spearman's
%                                       rank correlation between
%                                       observations (treated as sequences
%                                       of values).
%                                   'hamming'     
%                                     - Hamming distance, percentage of
%                                       coordinates that differ.
%                                   'jaccard'     
%                                     - One minus the Jaccard coefficient,
%                                       the percentage of nonzero
%                                       coordinates that differ.
%                                   function     
%                                     - A distance function specified using
%                                       @ (for example @DISTFUN). A
%                                       distance function must be of the
%                                       form
%
%                                       function D2 = DISTFUN(ZI, ZJ),
%
%                                       taking as arguments a 1-by-N vector
%                                       ZI containing a single row of X or
%                                       Y, an M2-by-N matrix ZJ containing
%                                       multiple rows of X or Y, and
%                                       returning an M2-by-1 vector of
%                                       distances D2, whose Jth element is
%                                       the distance between the
%                                       observations ZI and ZJ(J,:).
%        'Exponent'             - A positive scalar indicating the
%                                   exponent of Minkowski distance. This
%                                   argument is only valid when
%                                  'Distance' is 'minkowski'. Default: 2.
%        'Cov'                  - A positive definite matrix indicating the
%                                 covariance matrix when computing the
%                                 Mahalanobis distance. This argument is
%                                 only valid when 'Distance' is
%                                 'mahalanobis'. Default is NANCOV(X).
%       'Scale'                 - A vector S containing non-negative
%                                 values, with length equal to the number
%                                 of columns in X. Each  coordinate
%                                 difference between X and a query point is
%                                 divided by the corresponding element of
%                                 S. This argument is only valid when
%                                 'Distance' is 'seuclidean'. Default is
%                                 NANSTD(X).
%       'BucketSize'            - The maximum number of data points in the
%                                 leaf node of the kd-tree (default is 50).
%                                 This argument is only meaningful when
%                                 'NSMethod' is 'kdtree'.
%
%   See also ClassificationKNN, fitensemble, ClassificationKNN.fit.

            classreg.learning.FitTemplate.catchType(varargin{:});
            temp = classreg.learning.FitTemplate.make('KNN','type','classification',varargin{:});
         end
    end
    
    
    methods (Hidden)
        %Make this method hidden because it won't be documented
        function cmp = compact(this)
        %COMPACT Compact discriminant analysis.
        %   CMP=COMPACT(KNN) returns an object of class
        %   ClassificationKNN 
               
            cmp = this; %just return the class of full class
        end
        
    end
end %classdef

function distWgt =wgtFunc(dist,e)
% Compute scoring weights for knn given observation weights and distances.
% Both are N-by-K matrices. This is for quadratic inverse distance
% weights(1./dist.^e).

%take the factor of minDist to avoid overflow when summing the weights in the future
% the distance values has to be non-negative values.
minDist = min(dist,[],2); 
distNormalized = bsxfun(@rdivide,dist,minDist);
distNormalized(dist==0)=1; % to deal with zero distance values 
distWgt = 1./(distNormalized.^e);

end 


