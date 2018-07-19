function [coeff, score, latent, tsquared, explained] = pca(x,varargin)
%PCA Principal Component Analysis (PCA) on raw data.
%   COEFF = PCA(X) returns the principal component coefficients for the N
%   by P data matrix X. Rows of X correspond to observations and columns to
%   variables. Each column of COEFF contains coefficients for one principal
%   component. The columns are in descending order in terms of component
%   variance (LATENT). PCA, by default, centers the data and uses the
%   singular value decomposition algorithm. For the non-default options,
%   use the name/value pair arguments.
%   
%   [COEFF, SCORE] = PCA(X) returns the principal component score, which is
%   the representation of X in the principal component space. Rows of SCORE
%   correspond to observations, columns to components. The centered data
%   can be reconstructed by SCORE*COEFF'.
%
%   [COEFF, SCORE, LATENT] = PCA(X) returns the principal component
%   variances, i.e., the eigenvalues of the covariance matrix of X, in
%   LATENT.
%
%   [COEFF, SCORE, LATENT, TSQUARED] = PCA(X) returns Hotelling's T-squared
%   statistic for each observation in X. PCA uses all principal components
%   to compute the TSQUARED (computes in the full space) even when fewer
%   components are requested (see the 'NumComponents' option below). For
%   TSQUARED in the reduced space, use MAHAL(SCORE,SCORE).
%
%   [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = PCA(X) returns a vector
%   containing the percentage of the total variance explained by each
%   principal component.
%
%   [...] = PCA(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional
%   parameter name/value pairs to control the computation and handling of
%   special data types. Parameters are:
%   
%     'Algorithm' - Algorithm that PCA uses to perform the principal
%                   component analysis. Choices are:
%         'svd'   - Singular Value Decomposition of X (the default).
%         'eig'   - Eigenvalue Decomposition of the covariance matrix. It
%                   is faster than SVD when N is greater than P, but less
%                   accurate because the condition number of the covariance
%                   is the square of the condition number of X.
%
%     'Centered' - Indicator for centering the columns of X. Choices are: 
%         true   - The default. PCA centers X by subtracting off column
%                  means before computing SVD or EIG. If X contains NaN
%                  missing values, NANMEAN is used to find the mean with
%                  any data available.
%         false  - PCA does not center the data. In this case, the original
%                  data X can be reconstructed by X = SCORE*COEFF'. 
%
%     'Economy'  - Indicator for economy size output, when D the degrees of
%                  freedom is smaller than P. D, is equal to M-1, if data
%                  is centered and M otherwise. M is the number of rows
%                  without any NaNs if you use 'Rows', 'complete'; or the
%                  number of rows without any NaNs in the column pair that
%                  has the maximum number of rows without NaNs if you use
%                  'Rows', 'pairwise'. When D < P, SCORE(:,D+1:P) and
%                  LATENT(D+1:P) are necessarily zero, and the columns of
%                  COEFF(:,D+1:P) define directions that are orthogonal to
%                  X. Choices are:
%         true   - This is the default. PCA returns only the first D
%                  elements of LATENT and the corresponding columns of
%                  COEFF and SCORE. This can be significantly faster when P
%                  is much larger than D.
%         false  - PCA returns all elements of LATENT. Columns of COEFF and
%                  SCORE corresponding to zero elements in LATENT are
%                  zeros.
%
%     'NumComponents' - The number of components desired, specified as a
%                  scalar integer K satisfying 0 < K <= P. When specified,
%                  PCA returns the first K columns of COEFF and SCORE.
%
%     'Rows'     - Action to take when the data matrix X contains NaN
%                  values. Choices are:
%         'complete' - The default action. Observations with NaN values
%                      are removed before calculation. Rows of NaNs are
%                      inserted back into SCORE at the corresponding
%                      location.
%         'pairwise' - This option only applies when 'eig' method is used.
%                      If specified, PCA switches 'Algorithm' to 'eig'. The
%                      (I,J) element of the covariance matrix is computed
%                      using rows with no NaN values in columns I or J of
%                      X. Please note that the resulting covariance matrix
%                      may not be positive definite. In that case, PCA
%                      terminates with an error message.
%         'all'      - X is expected to have no missing values. All data
%                      are used, and execution will be terminated if NaN is
%                      found.
%                     
%     'Weights'  - Observation weights, a vector of length N containing all
%                  positive elements.
%
%     'VariableWeights' - Variable weights. Choices are:
%          - a vector of length P containing all positive elements.
%          - the string 'variance'. The variable weights are the inverse of
%            sample variance. If 'Centered' is set true at the same time,
%            the data matrix X is centered and standardized. In this case,
%            PCA returns the principal components based on the correlation
%            matrix.
%
%   Example:
%       load hald;
%       [coeff, score, latent, tsquared, explained] = pca(ingredients);
%
%   See also PCACOV, PCARES, BIPLOT, BARTTEST, CANONCORR, FACTORAN,
%   ROTATEFACTORS.
%
%   References:
%       [1] Jolliffe, I.T. Principal Component Analysis, 2nd ed., Springer,
%           2002.
%       [2] Krzanowski, W.J., Principles of Multivariate Analysis, Oxford
%           University Press, 1988.
%       [3] Seber, G.A.F., Multivariate Observations, Wiley, 1984. 
%       [4] Jackson, J.E., A User's Guide to Principal Components, Wiley, 
%           1988.

%   Copyright 2012 The MathWorks, Inc.
%   $Revision: 1.1.6.5.2.1 $ $Date: 2012/07/12 14:04:32 $


[n, p] = size(x);

% Parse arguments and check if parameter/value pairs are valid 
paramNames = {'Algorithm','Centered','Economy','NumComponents','Rows',...
    'Weights','VariableWeights'};
defaults   = {'svd',       true,      true,    p,           'complete',...
    ones(1,n) ,ones(1,p)};

[vAlgorithm, vCentered, vEconomy, vNumComponents, vRows,vWeights,...
    vVariableWeights,setFlag]...
    = internal.stats.parseArgs(paramNames, defaults, varargin{:});
% Validate String value for  Algorithm value
AlgorithmNames = {'svd','eig'};
vAlgorithm = internal.stats.getParamVal(vAlgorithm,AlgorithmNames,...
    '''Algorithm''');
% Validate boolean value for 'Centered' option
vCentered = internal.stats.parseOnOff(vCentered,'''Centered''');
% Validate boolean value for 'Economy' option
vEconomy = internal.stats.parseOnOff(vEconomy,'''Economy''');
% Validate the number of components option 'NumComponents'
if ~isempty(x) && ~internal.stats.isScalarInt(vNumComponents,1,p)
    error(message('stats:pca:WrongNumComponents',p));
end
% Validate value for 'Rows' option
RowsNames = {'complete','pairwise','all'};
vRows = internal.stats.getParamVal(vRows,RowsNames,'''Rows''');
if strcmp(vRows,'pairwise')  
    if setFlag.Algorithm && strcmp(vAlgorithm,'svd')    
    warning(message('stats:pca:NoPairwiseSVD'));
    end
    vAlgorithm = 'eig';
end
% Validate Weights Vectors
if isvector(vWeights) && isequal(numel(vWeights),n)
    vWeights = reshape(vWeights,1,n); % make sure it is a row vector
else
    error(message('stats:pca:WrongObsWeights', n));
end

if internal.stats.isString(vVariableWeights)
    WeightsNames = {'variance'};
    internal.stats.getParamVal(vVariableWeights,WeightsNames,...
        '''VariableWeights''');
    vVariableWeights = 1./classreg.learning.internal.wnanvar(x,vWeights,1);
elseif isnumeric(vVariableWeights) && isvector(vVariableWeights)...
        && (isequal(numel(vVariableWeights), p))
    vVariableWeights = reshape(vVariableWeights,1,p);
else
    error(message('stats:pca:WrongVarWeights', p));
end

if any(vWeights <= 0) || any(vVariableWeights <= 0)
    error(message('stats:pca:NonPositiveWeights'));
end
% end of checking input arguments


% Handling special empty case
if isempty(x)
    pOrZero = ~vEconomy * p;
    coeff = zeros(p, pOrZero); 
    coeff(1:p+1:end) = 1;
    score = zeros(n, pOrZero);
    latent = zeros(pOrZero, 1);
    tsquared = zeros(n, 1);
    return
end

nanIdx = isnan(x);
numNaN = sum(nanIdx, 2); % number of NaNs in each row
wasNaN = any(numNaN,2); % Rows that contain NaN

if strcmp(vRows,'all') 
    if any(wasNaN)
        error(message('stats:pca:RowsAll'));
    else
        vRows = 'complete';
    end
end

if strcmp(vRows,'complete')
    % Degrees of freedom (DOF) is n-1 if centered and n if not centered,
    % where n is the numer of rows without any NaN element.
    DOF = max(0,n-vCentered-sum(wasNaN));
elseif strcmp(vRows,'pairwise') 
    % DOF is the maximum number of element pairs without NaNs
    notNaN = double(~nanIdx);
    nanC = notNaN'*notNaN;
    nanC = nanC.*(~eye(p));
    DOF = max(nanC(:));
    DOF = DOF-vCentered;
end

% Center the data if 'Centered' is true.
if vCentered
    x = bsxfun(@minus,x,classreg.learning.internal.wnanmean(x, vWeights));
end

% Compute by EIG if no weights are given
if strcmp(vAlgorithm, 'eig')
    % Use EIG to compute.
    [coeff, eigValues] = localEIG(x, vCentered, vRows, vWeights,...
        vVariableWeights);
    
    % When 'Economy' value is true, nothing corresponding to zero
    % eigenvalues should be returned.
    if (DOF<p)
        if vEconomy
            coeff(:, DOF+1:p) = [];
            eigValues(DOF+1:p, :) = [];
        else % make sure they are zeros.
            eigValues(DOF+1:p, :) = 0;
        end
    end
    
    % Check if eigvalues are all postive
    if any(eigValues<0)
        error(message('stats:pca:CovNotPositiveSemiDefinite'));
    end
    
    if nargout > 1
        score = x/coeff';
        latent = eigValues; % Output Eigenvalues
        if nargout > 3
            tsquared = localTSquared(score, latent, n, p);
        end
    end
    
else % Use SVD to compute
    [U,sigma, coeff, wasNaN] = localSVD(x, n,...
        vEconomy, vWeights, vVariableWeights);
    if nargout > 1
        score =  bsxfun(@times,U,sigma');
        latent = sigma.^2./DOF;
        if nargout > 3
            tsquared = localTSquared(score,latent,DOF,p);
        end
        %Insert NaNs back
        if any(wasNaN)
            score = statinsertnan(wasNaN, score);
            if nargout >3
                tsquared = statinsertnan(wasNaN,tsquared);
            end
        end
    end
    
    if DOF < p
        % When 'Economy' value is true, nothing corresponding to zero
        % eigenvalues should be returned.
        if vEconomy
            coeff(:, DOF+1:end) = [];
            if nargout > 1
                score(:, DOF+1:end)=[];
                latent(DOF+1:end, :)=[];
            end
        elseif nargout > 1
        % otherwise, eigenvalues and corresponding outputs need to pad
        % zeros because svd(x,0) does not return columns of U corresponding
        % to components of (DOF+1):p.
            score(:, DOF+1:p) = 0;
            latent(DOF+1:p, 1) = 0;
        end
    end
end % end of EIG or SVD computation

% Calcuate the percentage of the total variance explained by each principal
% component.
if nargout > 4
    explained = 100*latent/sum(latent);
end

% Output only the first k principal components
if (vNumComponents<DOF)
    coeff(:, vNumComponents+1:end) = [];
    if nargout > 1
    score(:, vNumComponents+1:end) = [];
    end
end


% Enforce a sign convention on the coefficients -- the largest element in
% each column will have a positive sign.
[~,maxind] = max(abs(coeff), [], 1);
[d1, d2] = size(coeff);
colsign = sign(coeff(maxind + (0:d1:(d2-1)*d1)));
coeff = bsxfun(@times, coeff, colsign);
if nargout > 1
    score = bsxfun(@times, score, colsign); % scores = score
end

end % End of main function


%----------------Subfucntions--------------------------------------------

function [coeff, eigValues]=localEIG(x,vCentered, vRows,vWeights,...
        vVariableWeights)
% Compute by EIG. vRows are the options of handing NaN when compute
% covariance matrix

% Apply observation and variable weights
OmegaSqrt = sqrt(vWeights);
PhiSqrt = sqrt(vVariableWeights);
x = bsxfun(@times, x, OmegaSqrt');
x = bsxfun(@times, x, PhiSqrt);

xCov = ncnancov(x, vRows, vCentered);

[coeff, eigValueDiag] = eig(xCov);
[eigValues, idx] = sort(diag(eigValueDiag), 'descend');
coeff = coeff(:, idx);

coeff = bsxfun(@times, coeff,1./PhiSqrt');
end


function [U,sigma, coeff, wasNaN] = localSVD(x, n,...,
    vEconomy, vWeights, vVariableWeights)
% Compute by SVD. Weights are supplied by vWeights and vVariableWeights.

% Remove NaNs missing data and record location
[~,wasNaN,x] = statremovenan(x);
if n==1  % special case because statremovenan treats all vectors as columns
    wasNaN = wasNaN';
    x = x';
end

% Apply observation and variable weights
vWeights(wasNaN) = [];
OmegaSqrt = sqrt(vWeights);
PhiSqrt = sqrt(vVariableWeights);
x = bsxfun(@times, x, OmegaSqrt');
x = bsxfun(@times, x, PhiSqrt);
    
if vEconomy
    [U,sigma,coeff] = svd(x,'econ');
else
    [U,sigma, coeff] = svd(x, 0);
end

U = bsxfun(@times, U, 1./OmegaSqrt');
coeff = bsxfun(@times, coeff,1./PhiSqrt');

if n == 1     % sigma might have only 1 row
    sigma = sigma(1);
else
    sigma = diag(sigma);
end
end

function tsquared = localTSquared(score, latent,DOF, p)
% Subfunction to calulate the Hotelling's T-squared statistic. It is the
% sum of squares of the standardized scores, i.e., Mahalanobis distances.
% When X appears to have column rank < r, ignore components that are
% orthogonal to the data.

if isempty(score)
    tsquared = score;
    return;
end

r = min(DOF,p); % Max possible rank of x;
if DOF > 1
    q = sum(latent > max(DOF,p)*eps(latent(1)));
    if q < r
        warning(message('stats:pca:ColRankDefX', q)); 
    end
else
    q = 0;
end
standScores = bsxfun(@times, score(:,1:q), 1./sqrt(latent(1:q,:))');
tsquared = sum(standScores.^2, 2);
end

function c = ncnancov(x,Rows,centered)
%   C = NCNANCOV(X) returns X'*X/N, where N is the number of observations
%   after removing rows missing values. 
%
%   C = NCNANCOV(...,'pairwise') computes C(I,J) using rows with no NaN
%   values in columns I or J.  The result may not be a positive definite
%   matrix. C = NCNANCOV(...,'complete') is the default, and it omits rows
%   with any NaN values, even if they are not in column I or J.
%
%   C = NCNANCOV(...,true), C is normalized by N-1 if data X is already
%   centered. The default is false.

if nargin <2
    Rows = 'complete';
end

d = 0;
if nargin>2
    d =  d + centered;
end

idxnan = isnan(x);

[n, p] = size(x);


if ~any(any(idxnan))
    c = x'*x/(n-d);
elseif strcmp(Rows,'complete')
    nanrows = any(idxnan,2);
    xNotNaN = x((~nanrows),:);
    denom = max(0, (size(xNotNaN,1)-d) );
    c = xNotNaN'*xNotNaN/denom;
elseif strcmp(Rows,'pairwise')
    c = zeros(p,class(x));
    for i = 1:p
        for j = 1:i
            NonNaNRows = ~any(idxnan(:,[i, j]),2);
            denom = max(0,(sum(NonNaNRows)-d));
            c(i,j) = x(NonNaNRows,i)'*x(NonNaNRows,j)/denom;
        end
    end
    c = c + tril(c,-1)';
end
end