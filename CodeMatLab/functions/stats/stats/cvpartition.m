classdef cvpartition
    %CVPARTITION Create a cross-validation partition for data.
    %   An object of the CVPARTITION class defines a random partition on a
    %   set of data of a specified size.  This partition can be used to
    %   define test and training sets for validating a statistical model
    %   using cross-validation.
    %
    %   C = CVPARTITION(N,'Kfold',K) creates a CVPARTITION object C
    %   defining a random partition for K-fold cross-validation on N
    %   observations. The partition divides N observations into K disjoint
    %   subsamples (folds), chosen randomly but with roughly equal size.
    %   The default value of K is 10.
    %
    %   C = CVPARTITION(GROUP,'Kfold',K) creates a CVPARTITION object C
    %   defining a random partition for a stratified K-fold
    %   cross-validation. GROUP is a vector indicating the class
    %   information for each observation. GROUP can be a categorical
    %   variable, a numeric vector, a string array, or a cell array of
    %   strings. Each subsample has roughly equal size and roughly the same
    %   class proportions as in GROUP. CVPARTITION treats NaNs or empty
    %   strings in GROUP as missing values.
    %
    %   C = CVPARTITION(N,'Holdout',P) creates a CVPARTITION object C
    %   defining a random partition for holdout validation on N
    %   observations. This partition divides N observations into a training
    %   set and a test (holdout) set. P must be a scalar. When 0<P<1,
    %   CVPARTITION randomly selects approximately P*N observations for the
    %   test set. When P is an integer, it randomly selects P observations
    %   for the test set. The default value of P is 1/10.
    %
    %   C = CVPARTITION(GROUP,'Holdout',P) randomly partitions observations
    %   into a training set and a test set with stratification, using the
    %   class information in GROUP, i.e., both training and test sets have
    %   roughly the same class proportions as in GROUP.
    %
    %   C = CVPARTITION(N,'Leaveout') creates an object C defining a random
    %   partition for Leave-one-out cross-validation on N observations.
    %   Leave-one-out is a special case of K-fold in which the number of
    %   folds is equal to the number of observations N.
    %
    %   C = CVPARTITION(N,'Resubstitution') creates a CVPARTITION object C
    %   which does not partition the data. Both the training set and the
    %   test set contain all of the original N observations.
    %
    %   C has the following properties:
    %
    %      Type         The type of validation partition. It is 'kfold',
    %                   'holdout', 'leaveout' or 'resubstitution'.
    %      N            The number of observations (including observations
    %                   with missing GROUP values if GROUP is provided).
    %      NumTestSets  The number of test sets. Its value is the
    %                   number of folds in K-fold and Leave-one-out, 1 in
    %                   Holdout and Resubstitution.
    %      TrainSize    The size of each training set. It is a vector in
    %                   K-fold and Leave-one-out, a scalar in Holdout and
    %                   Resubstitution.
    %      TestSize     The size of each test set. It is a vector in K-fold
    %                   and Leave-one-out or a scalar in Holdout and
    %                   Resubstitution.
    %
    %   Example: Use 10-fold stratified cross-validation to compute the
    %   misclassification error for CLASSIFY on iris data.
    %
    %     load('fisheriris');
    %     CVO = cvpartition(species,'k',10);
    %     err = zeros(CVO.NumTestSets,1);
    %     for i = 1:CVO.NumTestSets
    %         trIdx = CVO.training(i);
    %         teIdx = CVO.test(i);
    %         ytest = classify(meas(teIdx,:),meas(trIdx,:),species(trIdx,:));
    %         err(i) = sum(~strcmp(ytest,species(teIdx)));
    %     end
    %     cvErr = sum(err)/sum(CVO.TestSize);
    %
    %   See also CVPARTITION/TEST, CVPARTITION/TRAINING, CVPARTITION/REPARTITION, CROSSVAL.
    
    %   Copyright 2007-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.9 $  $Date: 2011/10/15 02:06:17 $

    properties(GetAccess = 'public', SetAccess = 'protected')
        Type = '';
        N = [];
        NumTestSets = [];
        TrainSize = [];
        TestSize = [];
    end

    properties(GetAccess = 'protected', SetAccess = 'protected')
        indices = [];
        Group = [];
        holdoutT = [];
    end

    methods(Access = 'public')
        function cv = cvpartition(N,method,T,varargin)
            
            if isempty(varargin)
                s = RandStream.getGlobalStream;
                stdargin = nargin;
            else
                if length(varargin)>1
                    error(message('stats:cvpartition:cvpartTooManyArg'));
                else
                    stdargin = nargin-1;
                    s = varargin{1};
                    if ~isa(s,'RandStream')
                        error(message('stats:cvpartition:cvpartBadArg'));
                    end
                end
            end
            
            if stdargin < 2
                error(message('stats:cvpartition:TooFewInputs'));
            end

            if ischar(method) && size(method,1) == 1
                methodNames = {'kfold','holdout','leaveout','resubstitution'};
                j = find(strncmpi(method,methodNames,length(method)));
                if length(j) > 1
                    error(message('stats:cvpartition:AmbiguousMethod', method));
                elseif isempty(j)
                    error(message('stats:cvpartition:UnknownMethod'));
                end
            else
                error(message('stats:cvpartition:InvalidType'));
            end

            cv.Type = methodNames{j};

            if isscalar(N)
                if ~isnumeric(N) || N <= 1 || N ~= round(N) || ~isfinite(N)
                    error(message('stats:cvpartition:BadNX'));      
                end
                cv.N = N;
            else
                cv.Group = grp2idx(N);
                cv.N = length(cv.Group); % the number of observations including NaNs
                [~,wasnan,cv.Group] = statremovenan(cv.Group);
                hadNaNs = any(wasnan);
                if hadNaNs
                    warning(message('stats:cvpartition:MissingGroupsRemoved'));
                    if length (cv.Group) <= 1
                        error(message('stats:cvpartition:BadNGrp'));
                    end
                end
            end

            dftK = 10; % the default number of subsamples(folds) for Kfold
            P  = 1/10; % the default holdout ratio

            switch cv.Type
                case 'kfold'
                    if stdargin == 2 || isempty(T)
                        T = dftK;
                    elseif ~isscalar(T) || ~isnumeric(T) || T <= 1 || ...
                            T ~= round(T) || ~isfinite(T)
                        error(message('stats:cvpartition:BadK'));
                    end

                    if  isempty(cv.Group) && T > cv.N
                        warning(message('stats:cvpartition:KfoldGTN'));
                        T = cv.N;
                    elseif ~isempty(cv.Group) && T > length(cv.Group)
                        warning(message('stats:cvpartition:KfoldGTGN'));
                        T = length(cv.Group);
                    end

                    cv.NumTestSets = T; %set the number of fold
                    cv = cv.rerandom(s);

                case 'leaveout'
                    if stdargin > 2 && T ~= 1
                        error(message('stats:cvpartition:UnsupportedLeaveout'));
                    end

                    if isempty(cv.Group)
                        cv.NumTestSets = cv.N;
                    else
                        cv.NumTestSets = length(cv.Group);
                    end

                    [~,cv.indices] = sort(rand(s,cv.NumTestSets,1));

                    cv.TrainSize = (cv.NumTestSets-1) * ones(1,cv.NumTestSets);
                    cv.TestSize = ones(1,cv.NumTestSets);

                case 'resubstitution'
                    if stdargin > 2 
                        error(message('stats:cvpartition:UnsupportedCV'));
                    end

                    if isempty(cv.Group)
                        numObs = N;
                    else
                        numObs = length(cv.Group);
                    end

                    cv.indices = (1: numObs)';
                    cv.NumTestSets = 1;
                    cv.TrainSize =  numObs;
                    cv.TestSize =  numObs;

                case 'holdout'
                    if stdargin == 2 || isempty(T)
                        T = P;
                    elseif ~isscalar(T) || ~ isnumeric(T) || T <= 0 || ~isfinite(T)
                        error(message('stats:cvpartition:BadP'));
                    end

                    if T >= 1 %hold-T observations out
                        if T ~=round(T)
                            error(message('stats:cvpartition:BadP'));
                        end
                        if isempty(cv.Group)
                            if T >= cv.N
                                error(message('stats:cvpartition:PNotLTN'));
                            end
                        else
                            if T>= length(cv.Group)
                                error(message('stats:cvpartition:PNotLTGN'));
                            end
                        end
                    else
                        if (isempty(cv.Group) && floor(cv.N *T) == 0) ||...
                                (~isempty(cv.Group) && floor(length(cv.Group) * T) == 0)
                            error(message('stats:cvpartition:PTooSmall'));

                        end
                    end

                    cv.holdoutT = T;
                    cv = cv.rerandom(s);
                    cv.NumTestSets = 1;
            end

            %add NaNs back
            if ~isempty(cv.Group) && hadNaNs
                [cv.indices, cv.Group] =...
                    statinsertnan(wasnan, cv.indices, cv.Group);
            end
        end % cvpartition constructor

        function  cv = repartition(cv,varargin)
        %REPARTITION Rerandomize a cross-validation partition. 
        %   D = REPARTITION(C) creates a new random cross-validation partition D
        %   of the same type and size as C.  Use REPARTITION to perform multiple
        %   Monte-Carlo repetitions of cross-validation.
        %   D = REPARTITION(C,S) uses the RandStream object S as its
        %   random number generator.
        %   
        %   See also CVPARTITION.

            if isempty(varargin)
                s = RandStream.getGlobalStream;
            else
                if length(varargin)>1
                    error(message('stats:cvpartition:RepartTooManyArg'));
                else
                    s = varargin{1};
                    if ~isa(s,'RandStream')
                        error(message('stats:cvpartition:RepartBadArg'));
                    end
                end
            end
        
            if strcmp(cv.Type,'resubstitution')
                warning(message('stats:cvpartition:RepartNone'));
                return;
            end
            %remove NaNs from cv.Group
            if ~isempty(cv.Group)
                [~,wasnan,cv.Group] = statremovenan(cv.Group);
                hadNaNs = any(wasnan);
            end

            %  regenerate the data partition
            cv = cv.rerandom(s);

            %add NaNs back into cv.indices and cv.Group
            if ~isempty(cv.Group) && hadNaNs
                [cv.indices, cv.Group] =...
                    statinsertnan(wasnan, cv.indices, cv.Group);
            end
        end % repartition
       
        function trainIndices = training(cv,i)
        %TRAINING Training set for a cross-validation partition.
        %   TRIDX = TRAINING(C) returns a logical vector TRIDX that selects
        %   the observations in the training set for the hold-out
        %   cross-validation partition C.  C may also be a partition for
        %   resubstitution, in which case TRIDX is a logical vector that
        %   selects all observations.
        %
        %   TRIDX = TRAINING(C,I) returns a logical vector TRIDX that selects
        %   the observations in the I-th training set for a K-fold or
        %   leave-one-out cross-validation partition C.  In K-fold
        %   cross-validation, C divides a data set into K disjoint folds with
        %   roughly equal size.  The I-th training set consists of all
        %   observations not contained in the I-th fold.  In leave-one-out
        %   cross-validation, the I-th training set consists of the entire
        %   data set except the I-th observation.
        %
        %   See also CVPARTITION, CVPARTITION/TEST.
            switch cv.Type
                case {'kfold', 'leaveout'}
                    if nargin ~= 2
                        error(message('stats:cvpartition:WrongNumInputs'));
                    end
                    checkindex(i,cv.NumTestSets);
                    trainIndices = (cv.indices ~= i & ~isnan(cv.indices));

                case 'holdout'
                    if nargin == 2 && i~=1
                        error(message('stats:cvpartition:InvalidHOIndex'));
                    end
                    trainIndices = (cv.indices == 1);
                case 'resubstitution'
                    if nargin == 2 && i~= 1
                        error(message('stats:cvpartition:InvalidResubIndex'));
                    end
                    trainIndices = ~isnan(cv.indices);
            end
        end

        function testIndices = test(cv,i)
        %TEST Test set for a cross-validation partition.
        %   TEIDX = TEST(C) returns a logical vector TEIDX that selects the
        %   observations in the test set for the hold-out cross-validation
        %   partition C.  C may also be a partition for resubstitution, in
        %   which case TEIDX is a logical vector that selects all
        %   observations.
        %
        %   TEIDX = TEST(C,I) returns a logical vector TEIDX that selects the
        %   observations in the I-th test set for a K-fold or leave-one-out
        %   cross-validation partition C.  In K-fold cross-validation, C
        %   divides a data set into K disjoint folds with roughly equal size.
        %   The I-th test set consists of the I-th fold.  In leave-one-out
        %   cross-validation, the I-th test set consists of the I-th
        %   observation.
        %
        %   See also CVPARTITION, CVPARTITION/TRAINING.
            switch cv.Type
                case {'kfold','leaveout'}
                    if nargin ~= 2
                        error(message('stats:cvpartition:WrongNumInputs'));
                    end
                    checkindex(i,cv.NumTestSets);

                    testIndices = (cv.indices == i);

                case 'holdout'
                    if nargin == 2 && i~= 1
                        error(message('stats:cvpartition:InvalidHOIndex'));
                    end
                    testIndices = (cv.indices == 2);
                case 'resubstitution'
                    if nargin == 2 && i ~= 1
                        error(message('stats:cvpartition:InvalidResubIndex'));
                    end
                    testIndices = ~isnan(cv.indices);
            end
        end

        % Display methods
        function display(cv)
            isLoose = strcmp(get(0,'FormatSpacing'),'loose');

            objectname = inputname(1);
            if isempty(objectname)
                objectname = 'ans';
            end

            if (isLoose)
                fprintf('\n');
            end
            fprintf('%s = \n', objectname);
            disp(cv);
        end
        function disp(cv)
            isLoose = strcmp(get(0,'FormatSpacing'),'loose');

            if (isLoose)
                fprintf('\n');
            end
            switch cv.Type
                case 'kfold'
                    disp(getString(message('stats:cvpartition:KfoldCV')));
                case 'holdout'
                    disp(getString(message('stats:cvpartition:HoldoutCV')));
                case 'leaveout'
                    disp(getString(message('stats:cvpartition:LeaveoneoutCV')));
                case 'resubstitution'
                    disp (getString(message('stats:cvpartition:Resubstitution')));
            end
       
            fprintf('             N: %s\n',num2str(cv.N));
            fprintf('   NumTestSets: %s\n',num2str(cv.NumTestSets));
            Ndisp = 10;
            if cv.NumTestSets <= Ndisp
                 fprintf('     TrainSize: %s\n',num2str(cv.TrainSize ));
                 fprintf('      TestSize: %s\n',num2str(cv.TestSize ));

            else
                disp(['     TrainSize: ' num2str(cv.TrainSize(1:Ndisp)), ' ...']);
                disp(['      TestSize: ' num2str(cv.TestSize(1:Ndisp)), ' ...']);
            end
            %             end
        end
    end % public methods block


    methods(Access = 'private')
        %re-generate the data partition using the RandStream object s
        function cv = rerandom(cv,s)

            switch cv.Type
                case 'kfold'
                    if isempty(cv.Group)
                        if cv.NumTestSets == cv.N
                            %special case of K-fold -- loocv
                            [~,cv.indices] = sort(rand(s,cv.NumTestSets,1)); % randperm 
                            cv.TestSize = ones(1,cv.NumTestSets);
                        else
                            cv.indices = kfoldcv(cv.N,cv.NumTestSets,s);
                            cv.TestSize = accumarray(cv.indices,1)';
                        end
                    else
                        if cv.NumTestSets == length(cv.Group)
                            %special case of K-fold -- loocv
                            [~,cv.indices] = sort(rand(s,cv.NumTestSets,1)); % randperm
                            cv.TestSize = ones(1,cv.NumTestSets);
                        else
                            cv.indices = stra_kfoldcv(cv.Group,cv.NumTestSets,s);
                            cv.TestSize = accumarray(cv.indices,1)';
                        end
                    end

                    cv.TrainSize = size(cv.indices,1) - cv.TestSize;

                case 'holdout'
                    if cv.holdoutT >= 1
                        if isempty(cv.Group)
                            cv.indices = holdoutcv(cv.N,cv.holdoutT,s);
                            cv.TestSize = cv.holdoutT;
                            cv.TrainSize = cv.N-cv.TestSize;
                        else
                            cv.indices = stra_holdoutcv(cv.Group, cv.holdoutT/length(cv.Group), s);
                            cv.TestSize = sum(cv.indices == 2);
                            cv.TrainSize = sum(cv.indices == 1);
                        end
                    else %hold cv.holdoutT*N out
                        if isempty(cv.Group)
                            cv.TestSize = floor(cv.N * cv.holdoutT);
                            cv.TrainSize = cv.N-cv.TestSize;
                            cv.indices = holdoutcv(cv.N,cv.TestSize,s);
                        else
                            cv.indices = stra_holdoutcv(cv.Group,cv.holdoutT,s);
                            cv.TestSize = sum(cv.indices == 2);
                            cv.TrainSize = sum(cv.indices == 1);
                        end
                    end

                case 'leaveout'
                    [~,cv.indices] = sort(rand(s,cv.NumTestSets,1));
            end
        end
    end % private methods block
    
    methods(Hidden = true)
        function b = fieldnames(a)
            b = properties(a);
        end
        
        % Methods that we inherit, but do not want
        function a = fields(varargin),     throwUndefinedError(); end
        function a = ctranspose(varargin), throwUndefinedError(); end
        function a = transpose(varargin),  throwUndefinedError(); end
        function a = permute(varargin),    throwUndefinedError(); end
        function a = reshape(varargin),    throwUndefinedError(); end
        function a = cat(varargin),        throwNoCatError(); end
        function a = horzcat(varargin),    throwNoCatError(); end
        function a = vertcat(varargin),    throwNoCatError(); end
    end
    methods(Hidden = true, Static = true)
        function a = empty(varargin)
            error(message('stats:cvpartition:NoEmptyAllowed', upper( mfilename )));
        end
    end
   
end % classdef

function throwNoCatError()
error(message('stats:cvpartition:NoCatAllowed', upper( mfilename )));
end

function throwUndefinedError()
st = dbstack;
name = regexp(st(2).name,'\.','split');
error(message('stats:cvpartition:UndefinedFunction', name{ 2 }, mfilename));
end

%----------------------------------------------------
%stratified k-fold cross-validation
function cvid=stra_kfoldcv(group,nfold,s)
size_groups = accumarray(group(:),1);
if any(size_groups < nfold & size_groups > 0)
    warning(message('stats:cvpartition:KFoldMissingGrp'));
end
N = size(group,1);
cvid = 1 + mod((1:N)',nfold);
idrand = group + rand(s,N,1);
[~,idx] = sort(idrand);
cvid(idx) = cvid;
end

%----------------------------------------------------
%kfold cross-validation without stratification
function cvid = kfoldcv(N,nfold,s)
cvid = 1 + mod((1:N)',nfold);
[~,indices] = sort(rand(s,1,N)); % randperm
cvid = cvid(indices);
end

%-----------------------------------------------------
%holdout without stratification
function  idx= holdoutcv(N,num_test,s)
idx = 2*ones(N,1);
idx(1:N-num_test) = 1;
[~,indices] = sort(rand(s,1,N)); % randperm
idx = idx(indices);
end

%-----------------------------------------------------
%stratified holdout
function idx = stra_holdoutcv(group,test_ratio,s)
N = length(group);
size_groups = accumarray(group(:),1);
num_test = floor(size_groups * test_ratio);

test_diff = floor(N * test_ratio) - sum(num_test);
%add 1 for groups which are not in the test set
if any(num_test == 0)
    v=(num_test == 0);
    v(cumsum(v) > test_diff) = false;
    num_test(v) = num_test(v) + 1;
    test_diff = test_diff - sum(v);
end


if test_diff > 0
    ng= numel(size_groups);
    wasfull  =(num_test == size_groups);
    full_len = sum(wasfull);
    add = [ones(test_diff,1);zeros(ng - test_diff - full_len,1)];
    [~,indices] =  sort(rand(s,1,(ng-full_len))); % randperm
    add = add(indices);
    x = zeros(size(wasfull));
    x(~wasfull,:) = add;
    num_test = num_test + x;

end

if any(num_test == 0)
    warning(message('stats:cvpartition:HOTestZero'));
end

if any(num_test == size_groups)
    warning(message('stats:cvpartition:HOTrainingZero'));
end

idx = 2*ones(N,1);
for i = 1:numel(size_groups)
    g_idx = find(group == i);
    idx(g_idx(1:size_groups(i)-num_test(i))) = 1;
    [~,indices] = sort(rand(s,1,(size_groups(i)))); % randperm
    idx(g_idx) = idx(g_idx( indices ));
end

end

%-----------------------------
function checkindex(i,imax)
if ~(isnumeric(i) && isscalar(i) && i == round(i) && 1 <= i && i <= imax)
    error(message('stats:cvpartition:InvalidKFIndex'));
end
end

