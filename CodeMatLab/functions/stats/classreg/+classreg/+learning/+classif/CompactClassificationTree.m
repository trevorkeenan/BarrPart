classdef CompactClassificationTree < classreg.learning.classif.ClassificationModel
%CompactClassificationTree Decision tree for classification.
%   CompactClassificationTree is a decision tree with binary splits for
%   classification. It can predict response for new data.
%
%   CompactClassificationTree properties:
%       PredictorNames        - Names of predictors used for this tree.
%       CategoricalPredictors - Indices of categorical predictors.
%       ResponseName          - Name of the response variable.
%       ClassNames            - Names of classes in Y.
%       Cost                  - Misclassification costs.
%       Prior                 - Prior class probabilities.
%       ScoreTransform        - Transformation applied to predicted classification scores.
%       CatSplit              - Categorical splits for tree nodes.
%       Children              - Child nodes for tree nodes.
%       ClassCount            - Class counts for tree nodes.
%       ClassProb             - Class probabilities for tree nodes.
%       CutCategories         - Categories for splits on categorical predictors.
%       CutPoint              - Points for splits on continuous predictors.
%       CutType               - Split types (continuous or categorical) for tree nodes.
%       CutVar                - Split predictors for tree nodes.
%       IsBranch              - Branch (non-terminal) nodes.
%       NodeClass             - Majority class per tree node.
%       NodeErr               - Resubstitution error per tree node.
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
%   CompactClassificationTree methods:
%       edge                  - Classification edge.
%       loss                  - Classification loss.
%       margin                - Classification margins.
%       meanSurrVarAssoc      - Mean predictive measures of association between predictors based on surrogate splits.
%       predict               - Predicted response of this tree.
%       predictorImportance   - Importance of predictors for this tree.
%       view                  - View the tree structure.
%
%   See also ClassificationTree.

%   Copyright 2010-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.10 $  $Date: 2012/05/04 00:11:45 $
   
    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %CATSPLIT Categorical splits for tree nodes.
        %   The CatSplit property is an N-by-2 cell array, where N is the number of
        %   categorical splits in tree. Each row in CatSplit gives left and right
        %   values for a categorical split. For each branch node with categorical
        %   split J based on a predictor variable Z, the left child is chosen if Z
        %   is in CatSplit(J,1) and the right child is chosen if Z is in
        %   CatSplit(J,2). The splits are in the same order as nodes of the tree.
        %   Find the nodes for these splits by selecting 'categorical' cuts from
        %   top to bottom in the CutType property.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        CatSplit;
        
        %CHILDREN Child nodes for tree nodes.
        %   The Children property is an N-by-2 array containing the numbers of the
        %   child nodes for each node in tree, where n is the number of nodes. Leaf
        %   nodes have child node 0.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        Children;
        
        %CLASSCOUNT Class counts for tree nodes.
        %   The ClassCount property is an N-by-K array of class counts for the
        %   nodes in tree, where N is the number of nodes and K is the number of
        %   classes. For any node number I, the class counts ClassCount(I,:) are
        %   counts of observations (from the data used in fitting the tree) from
        %   each class satisfying the conditions for node I.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        ClassCount;
        
        %CLASSPROB Class probabilities for tree nodes.
        %   The ClassProb property is an N-by-K array of class probabilities for
        %   the nodes in tree, where N is the number of nodes and K is the number
        %   of classes. For any node number I, the class probabilities
        %   ClassProb(I,:) are the estimated probabilities for each class for a
        %   point satisfying the conditions for node I.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        ClassProb;
        
        %CUTCATEGORIES Categories for splits on categorical predictors.
        %   The CutCategories property is an N-by-2 cell array of the categories
        %   used at branches in the tree, where N is the number of nodes. For each
        %   branch node I based on a categorical predictor Z, the left child is
        %   chosen if Z is among the categories listed in CutCategories{I,1}, and
        %   the right child is chosen if Z is among those listed in
        %   CutCategories{I,2}. Both columns of CutCategories are empty for branch
        %   nodes based on continuous predictors and for leaf nodes.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        CutCategories;
        
        %CUTPOINT Points for splits on continuous predictors.
        %   The CutPoint property is an N-element vector of the values used as cut
        %   points in tree, where N is the number of nodes. For each branch node I
        %   based on a continuous predictor variable Z, the left child is chosen if
        %   Z < CutPoint(I) and the right child is chosen if Z >= CutPoint(I).
        %   CutPoint is NaN for branch nodes based on categorical predictors and
        %   for leaf nodes.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        CutPoint;
        
        %CUTTYPE Split types (continuous or categorical) for tree nodes.
        %   The CutType property is an N-element cell array indicating the type of
        %   cut at each node in the tree, where N is the number of nodes. For each
        %   node I, CutType{I} is:
        %       - 'continuous'  If the cut is defined in the form Z < V for a
        %                       variable Z and cut point V.
        %       - 'categorical' If the cut is defined by whether a variable Z
        %                       takes a value in a set of categories.
        %       - ''            If I is a leaf node.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        CutType;
        
        %CUTVAR Split predictors for tree nodes.
        %   The CutVar property is an N-element cell array of the names of the
        %   variables used for branching in each node in tree, where N is the
        %   number of nodes. For leaf nodes, CutVar contains an empty string.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        CutVar;
        
        %ISBRANCH Branch (non-terminal) nodes.
        %   The IsBranch property is an N-element logical vector that is true for
        %   each branch node and false for each leaf node of tree.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        IsBranch;
        
        %NODECLASS Majority class per tree node.
        %   The NodeClass property is an N-element cell array with the names of the
        %   most probable classes in each node of tree, where N is the number of
        %   nodes in the tree. Every element of this array is a string equal to one
        %   of the class names in ClassNames.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        NodeClass;
        
        %NODEERR Resubstitution error per tree node.
        %   The NodeErr property is an N-element vector of the errors of the nodes
        %   in the tree, where N is the number of nodes. NodeErr(I) is the
        %   misclassification probability for node I.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        NodeErr;
        
        %NODEPROB Tree node probability.
        %   The NodeProb property is an N-element vector of the probabilities of
        %   the nodes in tree, where N is the number of nodes. The probability of a
        %   node is computed as the proportion of observations from the original
        %   data that satisfy the conditions for the node. This proportion is
        %   adjusted for any prior probabilities assigned to each class.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        NodeProb;
        
        %NODESIZE Tree node size.
        %   The NodeSize property is an N-element vector of the sizes of the nodes
        %   in the tree, where N is the number of nodes. The size of a node is
        %   defined as the number of observations from the data used to create the
        %   tree that satisfy the conditions for the node.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        NodeSize;
        
        %NUMNODES Number of nodes in this tree.
        %   The NumNodes property is a numeric scalar showing the number of nodes
        %   in the tree.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        NumNodes;
        
        %PARENT Parent nodes for tree nodes.
        %   The Parent property is an N-element vector containing the number of the
        %   parent node for each node in tree, where N is the number of nodes. The
        %   parent of the root node is 0.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        Parent;
        
        %PRUNEALPHA Pruning values of alpha.
        %   The PruneAlpha property is a numeric vector with one element per
        %   pruning level. If the pruning level ranges from 0 to M, PruneAlpha has
        %   M+1 elements sorted in the ascending order: PruneAlpha(1) is for
        %   pruning level 0 (no pruning), PruneAlpha(2) is for pruning level 1 etc.
        %
        %   See also classreg.learning.classif.CompactClassificationTree,
        %   PruneList.
        PruneAlpha;
        
        %PRUNELIST Pruning sequence for this tree.
        %   The PruneList property is an N-element numeric vector with the pruning
        %   levels in each node of tree, where N is the number of nodes. The
        %   pruning levels range from 0 (no pruning) to M, where M is the distance
        %   between the deepest leaf and the root node.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        PruneList;
        
        %RISK Resubstitution risk per tree node.
        %   The Risk property is an N-element vector of the risk of the nodes in
        %   tree, where N is the number of nodes. If the tree was grown using an
        %   impurity measure (Gini index or deviance), the risk for each node is
        %   the measure of impurity for this node weighted by the node probability.
        %   Otherwise the risk for each node is the node error weighted by the node
        %   probability.
        %
        %   See also classreg.learning.classif.CompactClassificationTree,
        %   classreg.learning.classif.CompactClassificationTree/NodeErr,
        %   classreg.learning.classif.CompactClassificationTree/NodeProb.
        Risk;
        
        %SURRCUTCATEGORIES Categories for surrogate splits on categorical predictors.
        %   The SurrCutCategories property is an N-element cell array of the
        %   categories used for surrogate splits in the tree, where N is the number
        %   of nodes in tree. For each node I, SurrCutCategories{I} is a cell
        %   array. The length of SurrCutCategories{I} is equal to the number of
        %   surrogate predictors found at this node. Every element of
        %   SurrCutCategories{I} is either an empty string for a continuous
        %   surrogate predictor, or is a two-element cell array with categories for
        %   a categorical surrogate predictor. The first element of this
        %   two-element cell array lists categories assigned to the left child by
        %   this surrogate split and the second element of this two-element cell
        %   array lists categories assigned to the right child by this surrogate
        %   split. The order of the surrogate split variables at each node is
        %   matched to the order of variables in SurrCutVar. The optimal-split
        %   variable at this node does not appear. For non-branch (leaf) nodes,
        %   SurrCutCategories contains an empty cell.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        SurrCutCategories;
        
        %SURRCUTFLIP Signs of surrogate splits on continuous predictors.
        %   The SurrCutFlip property is an N-element cell array of the numeric cut
        %   assignments used for surrogate splits in the tree, where N is the
        %   number of nodes in tree. For each node I, SurrCutFlip{I} is a numeric
        %   vector. The length of SurrCutFlip{I} is equal to the number of
        %   surrogate predictors found at this node. Every element of
        %   SurrCutFlip{I} is either zero for a categorical surrogate predictor, or
        %   a numeric cut assignment for a continuous surrogate predictor. The
        %   numeric cut assignment can be either -1 or +1. For every surrogate
        %   split with a numeric cut C based on a continuous predictor variable Z,
        %   the left child is chosen if Z < C and the cut assignment for this
        %   surrogate split is +1, or if Z >= C and the cut assignment for this
        %   surrogate split is -1. Similarly, the right child is chosen if Z >= C
        %   and the cut assignment for this surrogate split is +1, or if Z < C and
        %   the cut assignment for this surrogate split is -1. The order of the
        %   surrogate split variables at each node is matched to the order of
        %   variables in SurrCutVar. The optimal-split variable at this node does
        %   not appear. For non-branch (leaf) nodes, SurrCutFlip contains an empty
        %   array.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        SurrCutFlip;
        
        %SURRCUTPOINT Points for surrogate splits on continuous predictors.
        %   The SurrCutPoint property is an N-element cell array of the numeric
        %   values used for surrogate splits in the tree, where N is the number of
        %   nodes in tree. For each node I, SurrCutPoint{I} is a numeric vector.
        %   The length of SurrCutPoint{I} is equal to the number of surrogate
        %   predictors found at this node. Every element of SurrCutPoint{I} is
        %   either NaN for a categorical surrogate predictor, or a numeric cut for
        %   a continuous surrogate predictor. For every surrogate split with a
        %   numeric cut C based on a continuous predictor variable Z, the left
        %   child is chosen if Z < C and SurrCutFlip for this surrogate split is
        %   -1. Similarly, the right child is chosen if Z >= C and SurrCutFlip for
        %   this surrogate split is +1, or if Z < C and SurrCutFlip for this
        %   surrogate split is -1. The order of the surrogate split variables at
        %   each node is matched to the order of variables returned by SurrCutVar.
        %   The optimal-split variable at this node does not appear. For non-branch
        %   (leaf) nodes, SurrCutPoint contains an empty cell.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        SurrCutPoint;
        
        %SURRCUTTYPE Surrogate split types (continuous or categorical) for tree nodes.
        %   The SurrCutType property is an N-element cell array indicating types of
        %   surrogate splits at each node in the tree, where N is the number of
        %   nodes in tree. For each node I, SurrCutType{I} is a cell array with the
        %   types of the surrogate split variables at this node. The variables are
        %   sorted by the predictive measure of association with the optimal
        %   predictor in the descending order, and only variables with the positive
        %   predictive measure are included. The order of the surrogate split
        %   variables at each node is matched to the order of variables in
        %   SurrCutVar. The optimal-split variable at this node does not appear.
        %   For non-branch (leaf) nodes, SurrCutType contains an empty cell. A
        %   surrogate split type can be either 'continuous' if the cut is defined
        %   in the form Z < V for a variable Z and cutpoint V or 'categorical' if
        %   the cut is defined by whether Z takes a value in a set of categories.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        SurrCutType;
        
        %SURRCUTVAR Surrogate split predictors for tree nodes.
        %   The SurrCutVar property is an N-element cell array of the names of the
        %   variables used for surrogate splits in each node in the tree, where N
        %   is the number of nodes in tree. Every element of SurrCutVar is a cell
        %   array with the names of the surrogate split variables at this node. The
        %   variables are sorted by the predictive measure of association with the
        %   optimal predictor in the descending order, and only variables with the
        %   positive predictive measure are included. The optimal-split variable at
        %   this node does not appear. For non-branch (leaf) nodes, SurrCutVar
        %   contains an empty cell.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        SurrCutVar;
        
        %SURRVARASSOC Predictive measures of association for surrogate splits for tree nodes.
        %   The SurrVarAssoc property is an N-element cell array of the predictive
        %   measures of association for surrogate splits in the tree, where N is
        %   the number of nodes in the tree. For each node I, SurrVarAssoc{I} is a
        %   numeric vector. The length of SurrVarAssoc{I} is equal to the number of
        %   surrogate predictors found at this node. Every element of
        %   SurrVarAssoc{I} gives the predictive measure of association between the
        %   optimal split and this surrogate split. The order of the surrogate
        %   split variables at each node is the order of variables in SurrCutVar.
        %   The optimal-split variable at this node does not appear. For non-branch
        %   (leaf) nodes, SurrVarAssoc contains an empty cell.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        SurrVarAssoc;
    end
    
    methods
        function a = get.CatSplit(this)
            a = catsplit(this.Impl.Tree);
        end
        
        function a = get.Children(this)
            a = children(this.Impl.Tree);
        end
        
        function a = get.ClassCount(this)
            a = classcount(this.Impl.Tree);
        end
        
        function a = get.ClassProb(this)
            a = classprob(this.Impl.Tree);
        end
        
        function a = get.CutCategories(this)
            a = cutcategories(this.Impl.Tree);
        end
        
        function a = get.CutPoint(this)
            a = cutpoint(this.Impl.Tree);
        end
        
        function a = get.CutType(this)
            a = cuttype(this.Impl.Tree);
        end
        
        function a = get.CutVar(this)
            a = cutvar(this.Impl.Tree);
        end
        
        function a = get.IsBranch(this)
            a = isbranch(this.Impl.Tree);
        end
        
        function a = get.NodeClass(this)
            a = nodeclass(this.Impl.Tree);
        end
        
        function a = get.NodeErr(this)
            a = nodeerr(this.Impl.Tree);
        end
        
        function a = get.NodeProb(this)
            a = nodeprob(this.Impl.Tree);
        end
        
        function a = get.NodeSize(this)
            a = nodesize(this.Impl.Tree);
        end
        
        function a = get.NumNodes(this)
            a = numnodes(this.Impl.Tree);
        end
        
        function a = get.Parent(this)
            a = parent(this.Impl.Tree);
        end
        
        function a = get.PruneAlpha(this)
            a = this.Impl.Tree.alpha;
        end
        
        function a = get.PruneList(this)
            a = prunelist(this.Impl.Tree);
        end
        
        function a = get.Risk(this)
            if isempty(this.Impl.Tree.impurity)
                a = risk(this.Impl.Tree,1:this.Impl.Tree.numnodes,'criterion','error');
            else
                a = risk(this.Impl.Tree,1:this.Impl.Tree.numnodes,'criterion','impurity');
            end
        end
        
        function a = get.SurrCutCategories(this)
            a = surrcutcategories(this.Impl.Tree);
        end
        
        function a = get.SurrCutFlip(this)
            a = surrcutflip(this.Impl.Tree);
        end
        
        function a = get.SurrCutPoint(this)
            a = surrcutpoint(this.Impl.Tree);
        end
        
        function a = get.SurrCutType(this)
            a = surrcuttype(this.Impl.Tree);
        end
        
        function a = get.SurrCutVar(this)
            a = surrcutvar(this.Impl.Tree);
        end
        
        function a = get.SurrVarAssoc(this)
            a = surrvarassoc(this.Impl.Tree);
        end        
    end
    
    methods(Access=protected)
        function this = CompactClassificationTree(dataSummary,classSummary,scoreTransform)
            this = this@classreg.learning.classif.ClassificationModel(...
                dataSummary,classSummary,scoreTransform);
        end
        
        function s = score(this,X,varargin)
            s = [];
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.ClassificationModel(this,s);
            s.CategoricalPredictors = this.CategoricalPredictors;
        end
    end
    
    methods
        function [labels,scores,node,cnum] = predict(this,X,varargin)
        %PREDICT Predict response of the tree.
        %   [LABEL,POSTERIOR,NODE,CNUM]=PREDICT(TREE,X) returns predicted class
        %   labels LABEL, posterior class probabilities POSTERIOR, node numbers
        %   NODE and class numbers CNUM for classification tree TREE and a matrix
        %   of predictors X. X must be a numeric matrix of size N-by-P, where P is
        %   the number of predictors used for training this model. Classification
        %   labels LABEL have the same type as Y used for training. Posterior class
        %   probabilities POSTERIOR are an N-by-K numeric matrix for N observations
        %   and K classes. NODE is a numeric vector of length N. CNUM is a numeric
        %   vector of length N with classes coded as integer numbers given by the
        %   order of classes in the ClassNames property. 
        %
        %   The predicted label is assigned to the class with the minimal
        %   misclassification cost. POSTERIOR*TREE.COST returns an N-by-K matrix of
        %   misclassification costs.
        %
        %   [LABEL,POSTERIOR,NODE,CNUM]=PREDICT(TREE,X,'PARAM1',val1,'PARAM2',val2,...)
        %   specifies optional parameter name/value pairs:
        %       'subtrees'    -  Vector SUBTREES of pruning levels, with 0
        %                        representing the full, unpruned tree. TREE must
        %                        include a pruning sequence as created by the
        %                        ClassificationTree.fit with 'prune' set to 'on' or
        %                        the PRUNE method. If SUBTREES has M elements and X
        %                        has N rows, then the output LABEL is an N-by-M
        %                        matrix, with the I-th column containing the fitted
        %                        values produced by the SUBTREES(I) subtree.
        %                        Similarly, POSTERIOR is an N-by-K-by-M matrix, and
        %                        NODE and CNUM are N-by-M matrices. SUBTREES must
        %                        be sorted in ascending order.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        
            % Get subtrees
            [subtrees,~,extraArgs] = ...
                internal.stats.parseArgs({'subtrees'},{0},varargin{:});
        
            % Check subtrees
            subtrees = processSubtrees(this.Impl,subtrees);

            % Get scores from the implementation
            if isempty(X)
                [labels,scores] = predictEmptyX(this,X);
                node = NaN(0,1);
                cnum = NaN(0,1);
                return;
            end
            [~,node] = eval(this.Impl.Tree,X,subtrees,extraArgs{:});            
            classnames = classreg.learning.internal.ClassLabel(this.Impl.Tree.classname);
            
            % Get the number of subtrees
            T = size(node,2);
            
            % Fill scores
            N = size(node,1);
            implscore = NaN(N,length(classnames),T);
            for t=1:T
                implscore(:,:,t) = classprob(this.Impl.Tree,node(:,t));
            end
            
            % Match the class name order
            K = length(this.ClassSummary.ClassNames);
            [~,pos] = ismember(classnames,this.ClassSummary.ClassNames);
            
            % Match scores to class names. Scores for trees represent
            % posterior class probabilities and by default are set to 0 for
            % missing classes.
            scores = zeros(N,K,T);
            scores(:,pos,:) = implscore;
 
            % Transform scores and find the most probable class
            prior = this.Prior;
            cost = this.Cost;
            scoreTransform = this.PrivScoreTransform;
            classnames = this.ClassNames;
            if ischar(classnames) && T>1
                classnames = cellstr(classnames);
            end
            labels = repmat(classnames(1,:),N,T);
            cnum = zeros(N,T);
            if T==1
                [labels,scores,~,cnum] = ...
                    this.LabelPredictor(classnames,prior,cost,scores,scoreTransform);
            else
                for t=1:T
                    [labels(:,t),scores(:,:,t),~,cnum(:,t)] = ...
                        this.LabelPredictor(classnames,prior,cost,scores(:,:,t),scoreTransform);
                end
            end
        end
        
        function [varargout] = loss(this,X,Y,varargin)
        %LOSS Classification error.
        %   L=LOSS(TREE,X,Y) returns classification cost for tree TREE computed
        %   using matrix of predictors X and true class labels Y. X must be a
        %   numeric matrix of size N-by-P, where P is the number of predictors used
        %   for training this model. Y must be of the same type as TREE.Y and have
        %   N elements.
        %
        %   L=LOSS(TREE,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'lossfun'   - Function handle for loss, or string representing a
        %                     built-in loss function. Available loss functions for
        %                     classification: 'binodeviance', 'classiferror',
        %                     'exponential', and 'mincost'. If you pass a function
        %                     handle FUN, LOSS calls it as shown below:
        %                          FUN(C,S,W,COST)
        %                     where C is an N-by-K logical matrix for N rows in X
        %                     and K classes in the ClassNames property, S is an
        %                     N-by-K numeric matrix, W is a
        %                     numeric vector with N elements, and COST is a K-by-K
        %                     numeric matrix. C has one true per row for the true
        %                     class. S is a matrix of posterior probabilities for
        %                     classes with one row per observation, similar to
        %                     POSTERIOR output from PREDICT. W is a vector of
        %                     observation weights. COST is a matrix of
        %                     misclassification costs. Default: 'mincost'
        %       'weights'   - Vector of observation weights. By default the weight
        %                     of every observation is set to 1. The length of this
        %                     vector must be equal to the number of rows in X.
        %
        %   [E,SE,NLEAF,BESTLEVEL]=LOSS(TREE,X,Y,'SUBTREES',SUBTREES) returns
        %   classification cost E, its standard error SE, number of leaves (terminal
        %   nodes) and the optimal pruning level for trees included in the pruning
        %   sequence SUBTREES. SUBTREES must be a vector with integer values between
        %   0 (full unpruned tree) and the maximal pruning level MAX(TREE.PruneList).
        %   E, SE and NLEAF are vectors of the same length as SUBTREES, and BESTLEVEL
        %   is a scalar. By default SUBTREES is set to 0. If you set SUBTREES to
        %   'all', LOSS uses the entire pruning sequence.
        %
        %   [...]=LOSS(TREE,X,Y,'SUBTREES',SUBTREES,'PARAM1',val1,'PARAM2',val2,...)
        %   specifies optional parameter name/value pairs:
        %      'treesize'   - If you choose 'se' (the default), LOSS returns
        %                     BESTLEVEL that corresponds to the smallest tree whose
        %                     cost is within one standard error of the minimum
        %                     cost. If you choose 'min', LOSS returns BESTLEVEL
        %                     that corresponds to the minimal cost tree.
        %
        %   See also classreg.learning.classif.CompactClassificationTree,
        %   classreg.learning.classif.CompactClassificationTree/predict.
            
            % Get input args
            N = size(X,1);
            args = {       'lossfun' 'subtrees' 'weights' 'treesize'};
            defs = {this.DefaultLoss          0 ones(N,1)       'se'};
            [funloss,subtrees,W,treesize] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check input args
            subtrees = processSubtrees(this.Impl,subtrees);
            funloss = classreg.learning.internal.lossCheck(funloss,'classification');
             
            % If user asks for the full tree, call the base method
            if isscalar(subtrees) && subtrees==0 && nargout<2
                [varargout{1:nargout}] = ...
                    loss@classreg.learning.classif.ClassificationModel(...
                    this,X,Y,'lossfun',funloss,'weights',W);
            else
                % Make sure the standard classification error is requested
                if ~strcmp(func2str(funloss),'classreg.learning.loss.mincost')
                    error(message('stats:classreg:learning:classif:CompactClassificationTree:loss:BadLossForSubtrees'));
                end
                [X,~,W,Y] = prepareDataForLoss(this,X,Y,W);
                [varargout{1:nargout}] = loss(this.Impl,X,cellstr(Y),'test',...
                    subtrees,treesize,'weights',W,...
                    'prior',this.Impl.Tree.prior,'cost',this.Impl.Tree.cost);
            end
        end
        
        function view(this,varargin)
        %VIEW View tree.
        %   VIEW(TREE) displays the decision tree.
        %
        %    VIEW(TREE,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %    parameter name/value pairs:
        %       'mode' - Either 'text' (default) or 'graph'. If 'text', VIEW
        %                displays TREE as text in the MATLAB command-line window.
        %                If 'graph', VIEW displays TREE as a figure in a new figure
        %                window.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
            
            view(this.Impl,varargin{:});
        end
                
        function [varargout] = meanSurrVarAssoc(this,varargin)
        %MEANSURRVARASSOC Mean predictive measure of association for surrogate splits in decision tree.
        %   MA=MEANSURRVARASSOC(TREE) returns a P-by-P matrix with predictive
        %   measures of association for P predictors. Element MA(I,J) is the
        %   predictive measure of association averaged over surrogate splits on
        %   predictor J for which predictor I is the optimal split predictor. This
        %   average is computed by summing positive values of the predictive
        %   measure of association over optimal splits on predictor I and surrogate
        %   splits on predictor J and dividing by the total number of optimal
        %   splits on predictor I, including splits for which the predictive
        %   measure of association between predictors I and J is negative.
        %
        %   MA=MEANSURRVARASSOC(T,N) takes an array N of node numbers and returns
        %   the predictive measure of association averaged over the specified
        %   nodes.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
            
            [varargout{1:nargout}] = meansurrvarassoc(this.Impl.Tree,varargin{:});
        end
        
        function imp = predictorImportance(this,varargin)
        %PREDICTORIMPORTANCE Estimates of predictor importance.
        %   IMP=PREDICTORIMPORTANCE(TREE) computes estimates of predictor importance
        %   for tree TREE by summing changes in the risk due to splits on every
        %   predictor and dividing the sum by the number of tree nodes. If the tree
        %   is grown without surrogate splits, this sum is taken over best splits
        %   found at each branch node. If the tree is grown with surrogate splits,
        %   this sum is taken over all splits at each branch node including
        %   surrogate splits. The returned vector IMP has one element for each
        %   input predictor in the data used to train this tree. At each node, the
        %   risk is estimated as node impurity if impurity was used to split nodes
        %   and node error otherwise. This risk is weighted by the node
        %   probability. Variable importance associated with this split is computed
        %   as the difference between the risk for the parent node and the total
        %   risk for the two children.
        %
        %   See also classreg.learning.classif.CompactClassificationTree.
        
            imp = varimportance(this.Impl.Tree,varargin{:});
        end
    end
        
end
