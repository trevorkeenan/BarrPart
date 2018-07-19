classdef TreeParams < classreg.learning.modelparams.ModelParams
%TreeParams Decision trees parameters.
%
%   TreeParams properties:
%       SplitCriterion    - 'gdi', 'twoing' or 'deviance'
%       MinParent         - Minimal size of parent node in tree
%       MinLeaf           - Minimal size of leaf node in tree
%       NVarToSample      - Number of predictors to select at random for decision split
%       MergeLeaves       - Flag for merging leaves after tree is grown
%       Prune             - Flag for computing the optimal pruning sequence
%       PruneCriterion    - 'error' (for classification and regression) or 'impurity' (for classification)
%       QEToler           - Tolerance on mean squared error per tree node (regression only)
%       Surrogate         - Flag for finding surrogate splits

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2011/05/09 01:23:45 $

    properties
        SplitCriterion = [];
        MinParent = [];
        MinLeaf = [];
        NVarToSample = [];
        MergeLeaves = [];
        Prune = [];
        PruneCriterion = [];
        QEToler = [];
        Surrogate = [];
    end

    methods(Access=protected)
        function this = TreeParams(type,splitcrit,minparent,...
                minleaf,nvartosample,mergeleaves,prune,prunecrit,qetoler,...
                surrogate)
            this = this@classreg.learning.modelparams.ModelParams('Tree',type);
            this.SplitCriterion = splitcrit;
            this.MinParent = minparent;
            this.MinLeaf = minleaf;
            this.NVarToSample = nvartosample;
            this.MergeLeaves = mergeleaves;
            this.Prune = prune;
            this.PruneCriterion = prunecrit;
            this.QEToler = qetoler;
            this.Surrogate = surrogate;
        end
    end

    methods(Static,Hidden)
        function [holder,extraArgs] = make(type,varargin)
            % Decode input args
            args = {'splitcriterion' 'minparent' 'minleaf' ...
                'nvartosample' 'mergeleaves' 'prune' 'prunecriterion' ...
                'qetoler' 'surrogate'};
            defs = {             []          []        [] ...
                            []            []      []               [] ...
                       []          []};
            [splitcrit,minparent,minleaf,...
                nvartosample,mergeleaves,prune,prunecrit,qetoler,...
                surrogate,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % No argument checking is performed here.
            % Normally, this would be the place to do it. Since
            % {Classification,Regression}Tree is just a wrapper to classregtree,
            % all parameter checks are carried out by classregtree.
            
            % Make argument holder
            holder = classreg.learning.modelparams.TreeParams(type,splitcrit,minparent,...
                minleaf,nvartosample,mergeleaves,prune,prunecrit,qetoler,...
                surrogate);
        end
    end

    methods(Hidden)
        function disp(this)
            fprintf('    SplitCriterion: %s\n',this.SplitCriterion);
            if ischar(this.MinParent)
                fprintf('         MinParent: %s\n',this.MinParent);
            else
                fprintf('         MinParent: %i\n',this.MinParent);
            end
            fprintf('           MinLeaf: %i\n',this.MinLeaf);
            if ischar(this.NVarToSample)
                fprintf('      NVarToSample: %s\n',this.NVarToSample);
            else
                fprintf('      NVarToSample: %i\n',this.NVarToSample);
            end
            fprintf('       MergeLeaves: %s\n',this.MergeLeaves);
            fprintf('             Prune: %s\n',this.Prune);
            fprintf('    PruneCriterion: %s\n',this.PruneCriterion);
            fprintf('         Surrogate: %s\n',this.Surrogate);

            isLoose = strcmp(get(0,'FormatSpacing'),'loose');
            if (isLoose)
                fprintf('\n');
            end
        end        
               
        function this = fillDefaultParams(this,X,Y,W,dataSummary,classSummary)
            if isempty(this.SplitCriterion)
                if strcmpi(this.Type,'classification')
                    this.SplitCriterion = 'gdi';
                else
                    this.SplitCriterion = 'mse';
                end
            end
            if isempty(this.MinParent)
                this.MinParent = 10;
            end
            if isempty(this.MinLeaf)
                this.MinLeaf = 1;
            end
            if isempty(this.NVarToSample)
                this.NVarToSample = 'all';
            end
            if isempty(this.MergeLeaves)
                this.MergeLeaves = 'on';
            end
            if isempty(this.Prune)
                this.Prune = 'on';
            end
            if isempty(this.PruneCriterion);
                this.PruneCriterion = 'error';
            end
            if strcmpi(this.Type,'regression') && isempty(this.QEToler)
                this.QEToler = 1e-6;
            end
            if isempty(this.Surrogate)
                this.Surrogate = 'off';
            end
        end
    end

end
