classdef RegressionTree < ...
        classreg.learning.regr.FullRegressionModel & classreg.learning.regr.CompactRegressionTree
%RegressionTree Decision tree for regression.
%   RegressionTree is a decision tree with binary splits for regression. It
%   can predict response for new data. It also stores data used for
%   training and can compute resubstitution predictions.
%
%   An object of this class cannot be created by calling the constructor.
%   Use RegressionTree.fit to create a RegressionTree object by fitting the
%   tree to training data.
%
%   This class is derived from CompactRegressionTree.
%
%   RegressionTree properties:
%       NObservations         - Number of observations.
%       X                     - Matrix of predictors used to train this tree.
%       Y                     - Observed response used to train this tree.
%       W                     - Weights of observations used to train this tree.
%       ModelParams           - Tree parameters.
%       PredictorNames        - Names of predictors used for this tree.
%       CategoricalPredictors - Indices of categorical predictors.
%       ResponseName          - Name of the response variable.
%       ResponseTransform     - Transformation applied to predicted regression response.
%       CatSplit              - Categorical splits for tree nodes.
%       Children              - Child nodes for tree nodes.
%       CutCategories         - Categories for splits on categorical predictors.
%       CutPoint              - Points for splits on continuous predictors.
%       CutType               - Split types (continuous or categorical) for tree nodes.
%       CutVar                - Split predictors for tree nodes.
%       IsBranch              - Branch (non-terminal) nodes.
%       NodeErr               - Resubstitution error per tree node.
%       NodeMean              - Mean response per tree node.
%       NodeProb              - Tree node probability.
%       NodeSize              - Tree node size.
%       NumNodes              - Number of nodes in this tree.
%       Parent                - Parent nodes for tree nodes.
%       PruneAlpha            - Pruning values of alpha.
%       PruneList             - Pruning sequence for this tree.
%       Risk                  - Resubstitution risk per tree node.
%       SurrCutCategories     - Categories for surrogate splits on categorical predictors.
%       SurrCutFlip           - Signs of surrogate splits on continuous predictors.
%       SurrCutPoint          - Points for surrogate splits on continuous predictors.
%       SurrCutType           - Surrogate split types (continuous or categorical) for tree nodes.
%       SurrCutVar            - Surrogate split predictors for tree nodes.
%       SurrVarAssoc          - Predictive measures of association for surrogate splits for tree nodes.
%
%   RegressionTree methods:
%       fit (static)          - Grow a regression tree.
%       template (static)     - Tree template for FITENSEMBLE.
%       compact               - Compact this tree.
%       crossval              - Cross-validate this tree.
%       cvLoss                - Regression loss by cross-validation.
%       loss                  - Regression loss.
%       meanSurrVarAssoc      - Mean predictive measures of association between predictors based on surrogate splits.
%       predict               - Predicted response of this tree.
%       predictorImportance   - Importance of predictors for this tree.
%       prune                 - Prune this tree.
%       resubLoss             - Resubstitution regression loss.
%       resubPredict          - Resubstitution predicted response.
%       view                  - View the tree structure.
%
%   Example: Grow a regression tree for car data.
%       load carsmall
%       t = RegressionTree.fit([Weight Horsepower],MPG,'predictornames',{'Weight' 'Horsepower'})
%       view(t)
%
%   See also classreg.learning.regr.CompactRegressionTree.
    
%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.10 $  $Date: 2012/01/27 21:05:19 $

    methods(Hidden)
        function this = RegressionTree(X,Y,W,modelParams,dataSummary,responseTransform)
            if nargin~=6 || ischar(W)
                error(message('stats:RegressionTree:RegressionTree:DoNotUseConstructor'));
            end
            this = this@classreg.learning.regr.FullRegressionModel(...
                X,Y,W,modelParams,dataSummary,responseTransform);
            this = this@classreg.learning.regr.CompactRegressionTree(...
                dataSummary,responseTransform);
            this = fitTree(this);
        end
    end

    methods(Static)
        function this = fit(X,Y,varargin)
        %FIT Fit a regression decision tree.
        %   TREE=REGRESSIONTREE.FIT(X,Y) returns a regression decision tree for
        %   predictors X and response Y.
        %
        %   X must be an N-by-P matrix of predictors with one row per observation
        %   and one column per predictor. Y must be a vector of floating-point
        %   numbers with N elements.
        %
        %   FIT grows the tree using MSE (mean squared error) as the splitting
        %   criterion.
        %
        %   TREE is a regression tree with binary splits. If you use one of the
        %   following five options, TREE is of class RegressionPartitionedModel:
        %   'crossval', 'kfold', 'holdout', 'leaveout' or 'cvpartition'. Otherwise,
        %   TREE is of class RegressionTree.
        %
        %   TREE=REGRESSIONTREE.FIT(X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %       'CategoricalPredictors' - List of categorical predictors. Pass
        %                                 'CategoricalPredictors' as one of:
        %                                   * A numeric vector with indices between
        %                                     1 and P, where P is the number of
        %                                     columns of X.
        %                                   * A logical vector of length P, where a
        %                                     true entry means that the
        %                                     corresponding column of X is a
        %                                     categorical variable.
        %                                   * 'all', meaning all predictors are
        %                                     categorical.
        %                                   * A cell array of strings, where each
        %                                     element in the array is the name of a
        %                                     predictor variable. The names must
        %                                     match entries in 'PredictorNames'
        %                                     values.
        %                                   * A character matrix, where each row of
        %                                     the matrix is a name of a predictor
        %                                     variable. The names must match
        %                                     entries in 'PredictorNames' values.
        %                                     Pad the names with extra blanks so
        %                                     each row of the character matrix has
        %                                     the same length.
        %                                 Default: []
        %       'crossval'              - If 'on', grows a cross-validated tree
        %                                 with 10 folds. You can use 'kfold',
        %                                 'holdout', 'leaveout' and 'cvpartition'
        %                                 parameters to override this
        %                                 cross-validation setting. You can only
        %                                 use one of these four options ('kfold',
        %                                 'holdout', 'leaveout' and 'cvpartition')
        %                                 at a time when creating a cross-validated
        %                                 tree. As an alternative, you can
        %                                 cross-validate later using CROSSVAL
        %                                 method for tree TREE. Default: 'off'
        %       'cvpartition'           - A partition created with CVPARTITION to
        %                                 use in cross-validated tree. 
        %       'holdout'               - Holdout validation uses the specified
        %                                 fraction of the data for test, and uses
        %                                 the rest of the data for training.
        %                                 Specify a numeric scalar between 0 and 1.
        %       'kfold'                 - Number of folds to use in cross-validated
        %                                 tree, a positive integer. Default: 10
        %       'leaveout'              - Use leave-one-out cross-validation by
        %                                 setting to 'on'. 
        %       'MergeLeaves'           - When 'on', RegressionTree.fit merges
        %                                 leaves that originate from the same
        %                                 parent node, and that give a sum of risk
        %                                 values greater or equal to the risk
        %                                 associated with the parent node. When
        %                                 'off', RegressionTree.fit does not
        %                                 merge leaves. Default: 'on'
        %       'MinLeaf'               - Each leaf has at least 'MinLeaf'
        %                                 observations per tree leaf. If you supply
        %                                 both 'MinParent' and 'MinLeaf',
        %                                 RegressionTree.fit uses the setting
        %                                 that gives larger leaves: MinParent =
        %                                 max(MinParent,2*MinLeaf). Default: 1
        %       'MinParent'             - Each splitting node in the tree has at
        %                                 least 'MinParent' observations. If you
        %                                 supply both 'MinParent' and 'MinLeaf',
        %                                 RegressionTree.fit uses the setting
        %                                 that gives larger leaves: MinParent =
        %                                 max(MinParent,2*MinLeaf). Default: 10
        %       'NVarToSample'          - Number of predictors to select at random
        %                                 for each split. Can be a positive integer
        %                                 or 'all'; 'all' means use all available
        %                                 predictors. Default: 'all'
        %       'PredictorNames'        - A cell array of names for the predictor
        %                                 variables, in the order in which they
        %                                 appear in X. Default: {'x1','x2',...}
        %       'prune'                 - When 'on', RegressionTree.fit grows
        %                                 the unpruned tree and computes the
        %                                 optimal sequence of pruned subtrees.
        %                                 When 'off' RegressionTree.fit grows
        %                                 the tree without pruning. Default: 'on'
        %       'ResponseName'          - Name of the response variable Y, a
        %                                 string. Default: 'Y'
        %       'surrogate'             - When 'on', RegressionTree.fit finds
        %                                 surrogate splits at each branch node.
        %                                 This setting can use much time and
        %                                 memory. You can use it to improve the
        %                                 tree accuracy for data with missing
        %                                 values. You can also use it to compute
        %                                 measures of predictive association
        %                                 between predictors. Default: 'off'
        %       'weights'               - Vector of observation weights, one weight
        %                                 per observation. RegressionTree.fit
        %                                 normalizes the weights to add up to one.
        %                                 Default: ones(size(X,1),1)
        %
        %   See also RegressionTree,
        %   classreg.learning.partition.RegressionPartitionedModel. 

            temp = RegressionTree.template(varargin{:});
            this = fit(temp,X,Y);
        end
        
        function temp = template(varargin)
        %TEMPLATE Create a regression template.
        %   T=REGRESSIONTREE.TEMPLATE() returns a tree template suitable for
        %   use in the FITENSEMBLE function.
        %
        %   T=REGRESSIONTREE.TEMPLATE('PARAM1',val1,'PARAM2',val2,...)
        %   specifies optional parameter name/value pairs:
        %       'MergeLeaves'           - When 'on', RegressionTree merges
        %                                 leaves that originate from the same
        %                                 parent node, and that give a sum of mean
        %                                 squared error (MSE) values greater or
        %                                 equal to the MSE associated with the
        %                                 parent node. When 'off', RegressionTree
        %                                 does not merge leaves. Default for
        %                                 ensemble fitting: 'off'
        %       'MinLeaf'               - Each leaf has at least 'MinLeaf'
        %                                 observations per tree leaf. If you supply
        %                                 both 'MinParent' and 'MinLeaf',
        %                                 RegressionTree uses the setting that
        %                                 gives larger leaves: MinParent =
        %                                 max(MinParent,2*MinLeaf). Default for
        %                                 ensemble fitting: 1 for boosting and 5
        %                                 for bagging.
        %       'MinParent'             - Each splitting node in the tree has at
        %                                 least 'MinParent' observations. If you
        %                                 supply both 'MinParent' and 'MinLeaf',
        %                                 RegressionTree uses the setting that
        %                                 gives larger leaves: MinParent =
        %                                 max(MinParent,2*MinLeaf). Default for
        %                                 ensemble fitting: number of training
        %                                 observations for boosting and ten for
        %                                 bagging.
        %       'NVarToSample'          - Number of predictors to select at random
        %                                 for each split. Can be a positive integer
        %                                 or 'all'; 'all' means use all available
        %                                 predictors. Default for ensemble fitting:
        %                                 'all' for boosting and one third of
        %                                 predictors for bagging.
        %       'prune'                 - When 'on', RegressionTree grows the
        %                                 unpruned tree and computes the optimal
        %                                 sequence of pruned subtrees.  When 'off'
        %                                 RegressionTree grows the tree without
        %                                 pruning. Default for ensemble fitting:
        %                                 'off'
        %       'PruneCriterion'        - String with the pruning criterion. The
        %                                 only criterion available for regression
        %                                 is 'mse' (mean squared error).
        %       'surrogate'             - When 'on', RegressionTree finds
        %                                 surrogate splits at each branch node.
        %                                 This setting can use much time and
        %                                 memory. You can use it to improve the
        %                                 tree accuracy for data with missing
        %                                 values. You can also use it to compute
        %                                 measures of predictive association
        %                                 between predictors. Default for ensemble
        %                                 fitting: 'off'
        %
        %   See also RegressionTree, fitensemble, RegressionTree.fit.            
            
            classreg.learning.FitTemplate.catchType(varargin{:});
            temp = classreg.learning.FitTemplate.make('Tree','type','regression',varargin{:});
        end
    end

    methods(Access=protected)
        function this = fitTree(this)
            if strcmpi(this.ModelParams.MinParent,'OneSplit')
                minparent = size(this.X,1);
            else
                minparent = this.ModelParams.MinParent;
            end
            tree = classregtree(this.X,this.Y,'weights',this.W,...
                'method','regression',...
                'names',this.DataSummary.PredictorNames,...
                'categorical',this.DataSummary.CategoricalPredictors,...
                'minparent',minparent,...
                'minleaf',this.ModelParams.MinLeaf,...
                'nvartosample',this.ModelParams.NVarToSample,...
                'mergeleaves',this.ModelParams.MergeLeaves,...
                'prune',this.ModelParams.Prune,...
                'qetoler',this.ModelParams.QEToler,...
                'surrogate',this.ModelParams.Surrogate,...
                'skipchecks',true);
            this.Impl = classreg.learning.impl.CompactTreeImpl(tree);
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.regr.CompactRegressionTree(this,s);
            s = propsForDisp@classreg.learning.regr.FullRegressionModel(this,s);
        end
    end
       
    methods
        function cmp = compact(this,varargin)
        %COMPACT Compact tree.
        %   CMP=COMPACT(TREE) returns an object of class CompactRegressionTree
        %   holding the structure of the trained tree. The compact object does not
        %   contain X and Y used for training.
        %
        %   See also RegressionTree,
        %   classreg.learning.regr.CompactRegressionTree.
        
            cmp = classreg.learning.regr.CompactRegressionTree(...
                this.DataSummary,this.PrivResponseTransform);
            cmp.Impl = this.Impl;
        end
 
        function this = prune(this,varargin)
        %PRUNE Produce a sequence of subtrees by pruning.
        %   T2=PRUNE(T1) computes the pruning sequence for decision tree T1 and
        %   returns a decision tree T2 with the pruning sequence filled in. Trees
        %   are pruned based on an optimal pruning scheme that first prunes
        %   branches giving less improvement in mean squared error. Subsequent
        %   calls to PRUNE can use this pruning sequence to reduce the tree by
        %   turning some branch nodes into leaf nodes, and removing the leaf nodes
        %   under the original branch.
        %
        %   T2=PRUNE(...,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs. You can use only one option at a time.
        %       'level'     - A numeric scalar between 0 (no pruning) and the
        %                     largest pruning level of this tree MAX(T1.PruneList).
        %                     PRUNE returns the tree pruned to this level.
        %       'nodes'     - A numeric vector with elements between 1 and
        %                     T1.NumNodes. Any T1 branch nodes listed in NODES
        %                     become leaf nodes in T2, unless their parent nodes
        %                     are also pruned.
        %       'alpha'     - A numeric positive scalar. PRUNE prunes T1 to the
        %                     specified value of the pruning cost.
        %
        %     T2=PRUNE(T1) returns the decision tree T2 that is the full, unpruned
        %     T1, but with optimal pruning information added.  This is useful only
        %     if you created T1 by pruning another tree, or by using the
        %     RegressionTree.fit with pruning set 'off'.  If you plan to prune
        %     a tree multiple times along the optimal pruning sequence, it is more
        %     efficient to create the optimal pruning sequence first.
        %
        %     Example:  Display full tree for car data, as well as the tree pruned
        %     to level 10.
        %        load carsmall;
        %        varnames = {'Weight' 'Horsepower'};
        %        t1 = RegressionTree.fit([Weight Horsepower],MPG,'predictornames',varnames)
        %        view(t1,'mode','graph');
        %        t2 = prune(t1,'level',10);
        %        view(t2,'mode','graph');
        %
        %   See also RegressionTree, PruneList, PruneAlpha.
            
            this.Impl = classreg.learning.impl.CompactTreeImpl(...
                prune(this.Impl.Tree,'criterion',...
                this.ModelParams.PruneCriterion,varargin{:}));
        end
        
        function [varargout] = resubPredict(this,varargin)
        %RESUBPREDICT Predict resubstitution response of the tree.
        %   [YFIT,NODE]=RESUBPREDICT(TREE) returns predicted response YFIT and node
        %   numbers NODE for regression tree TREE and training data TREE.X. YFIT is
        %   a vector of type double with NObservations elements. NODE is a numeric
        %   vector of length N.
        %
        %   [YFIT,NODE]=RESUBPREDICT(TREE,'PARAM1',val1,'PARAM2',val2,...)
        %   specifies optional parameter name/value pairs:
        %       'subtrees'    -  Vector SUBTREES of pruning levels, with 0
        %                        representing the full, unpruned tree. TREE must
        %                        include a pruning sequence as created either by
        %                        the RegressionTree.fit with 'prune' set to 'on',
        %                        or by the PRUNE method. If SUBTREES has M elements
        %                        and TREE.X has N rows, then the output YFIT is an
        %                        N-by-M matrix, with the I-th column containing the
        %                        fitted values produced by the SUBTREES(I) subtree.
        %                        Similarly, NODE is an N-by-M matrix. SUBTREES must
        %                        be sorted in ascending order.
        %
        %   See also RegressionTree, predict.
            
            [varargout{1:nargout}] = ...
                resubPredict@classreg.learning.regr.FullRegressionModel(this,varargin{:});
        end
        
        function [varargout] = resubLoss(this,varargin)
        %RESUBLOSS Regression error by resubstitution.
        %   L=RESUBLOSS(TREE) returns mean squared error for tree TREE
        %   computed for training data TREE.X and TREE.Y.
        %
        %   L=RESUBLOSS(TREE,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'lossfun'   - Function handle for loss, or string representing a
        %                     built-in loss function. Available loss functions for
        %                     regression: 'mse'. If you pass a function handle FUN,
        %                     LOSS calls it as shown below:
        %                          FUN(Y,Yfit,W)
        %                     where Y, Yfit and W are numeric vectors of length N.
        %                     Y is observed response, Yfit is predicted response,
        %                     and W is observation weights. Default: 'mse'
        %
        %   [E,SE,NLEAF,BESTLEVEL]=RESUBLOSS(TREE,'SUBTREES',SUBTREES) returns
        %   resubstitution error E, its standard error SE, number of leaves
        %   (terminal nodes) NLEAF and the optimal pruning level BESTLEVEL for
        %   trees included in the pruning sequence SUBTREES. SUBTREES must be a
        %   vector with integer values between 0 (full unpruned tree) and the
        %   maximal pruning level MAX(TREE.PruneList). E, SE and NLEAF are vectors
        %   of the same length as SUBTREES, and BESTLEVEL is a scalar. By default
        %   SUBTREES is set to 0. If you set SUBTREES to 'all', LOSS uses the
        %   entire pruning sequence.
        %
        %   [...]=RESUBLOSS(TREE,X,Y,'SUBTREES',SUBTREES,'PARAM1',val1,'PARAM2',val2,...)
        %   specifies optional parameter name/value pairs:
        %      'treesize'   - Either 'se' (the default) to choose the smallest
        %                     tree whose MSE is within one standard error of the
        %                     minimum MSE, or 'min' to choose the minimal MSE tree.
        %
        %   See also RegressionTree, loss.
        
            [varargout{1:nargout}] = ...
                resubLoss@classreg.learning.regr.FullRegressionModel(this,varargin{:});
        end
        
        function [varargout] = cvLoss(this,varargin)
        %CVLOSS Regression error by cross-validation.
        %   [E,SE,NLEAF,BESTLEVEL]=CVLOSS(TREE) returns cross-validated regression
        %   error E, its standard error SE, number of leaves (terminal nodes) NLEAF
        %   and the optimal pruning level BESTLEVEL for tree TREE.
        %
        %   [E,SE,NLEAF,BESTLEVEL]=CVLOSS(TREE,'PARAM1',val1,'PARAM2',val2,...)
        %   specifies optional parameter name/value pairs:
        %       'subtrees'    -  Vector SUBTREES of pruning levels, with 0
        %                        representing the full, unpruned tree. TREE must
        %                        include a pruning sequence as created either by
        %                        the RegressionTree.fit with 'prune' set to
        %                        'on', or by the PRUNE method. The returned E, SE
        %                        and NLEAF are vectors of the same length as
        %                        SUBTREES, and BESTLEVEL is a scalar. By default
        %                        SUBTREES is set to 0.  If you set SUBTREES to
        %                        'all', the entire pruning sequence will be used.
        %      'treesize'     -  If you choose 'se' (the default), CVLOSS returns
        %                        BESTLEVEL that corresponds to the smallest tree
        %                        whose MSE is within one standard error of the
        %                        minimum MSE. If you choose 'min', CVLOSS returns
        %                        BESTLEVEL that corresponds to the minimal MSE
        %                        tree. MSE is mean squared error.
        %      'kfold'        -  Number of cross-validation samples (default 10).
        %
        %   See also RegressionTree, loss.
            
            % Get input args
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            args = {'subtrees' 'kfold' 'treesize'};
            defs = {         0      10       'se'};
            [subtrees,kfold,treesize] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check input args
            subtrees = processSubtrees(this.Impl,subtrees);
            
            % Call classregtree/test
            [varargout{1:nargout}] = loss(this.Impl,this.X,this.Y,...
                'crossvalidate',subtrees,treesize,'nsamples',kfold,'weights',this.W);
        end
    end
end
