function [b,logL,H,stats] = coxphfit(X,y,varargin)
%COXPHFIT Cox proportional hazards regression.
%   B=COXPHFIT(X,Y) fits Cox's proportional hazards regression model to
%   the data in the vector Y, using the columns of the matrix X as predictors.
%   X and Y must have the same number of rows, and X should not contain a
%   column of ones.  The result B is a vector of coefficient estimates.
%   This model states that the hazard rate for the distribution of Y can be
%   written as h(t)*exp(X*B) where h(t) is a common baseline hazard function.
%
%   [...] = COXPHFIT(X,Y,'PARAM1',VALUE1,'PARAM2',VALUE2,...) specifies
%   additional parameter name/value pairs chosen from the following:
%
%      Name          Value
%      'baseline'    The X values at which the baseline hazard is to be
%                    computed.  Default is mean(X), so the hazard at X is
%                    h(t)*exp((X-mean(X))*B).  Enter 0 to compute the
%                    baseline relative to 0, so the hazard at X is
%                    h(t)*exp(X*B).
%      'censoring'   A boolean array of the same size as Y that is 1 for
%                    observations that are right-censored and 0 for
%                    observations that are observed exactly.  Default is
%                    all observations observed exactly.
%      'frequency'   An array of the same size as Y containing non-negative
%                    integer counts.  The jth element of this vector
%                    gives the number of times the jth element of Y and
%                    the jth row of X were observed.  Default is 1
%                    observation per row of X and Y.
%      'init'        A vector containing initial values for the estimated
%                    coefficients B.
%      'options'     A structure specifying control parameters for the
%                    iterative algorithm used to estimate B.  This argument
%                    can be created by a call to STATSET.  For parameter
%                    names and default values, type STATSET('coxphfit').
%
%   [B,LOGL,H,STATS]=COXPHFIT(...) returns additional results.  LOGL is the
%   log likelihood.  H is a two-column matrix containing y values in column
%   1 and the estimated baseline cumulative hazard evaluated at those
%   values in column 2.  STATS is a structure with the following fields:
%       'beta'     coefficient estimates (same as B output)
%       'se'       standard errors of coefficient estimates
%       'z'        z statistics for B (B divided by standard error)
%       'p'        p-values for B
%       'covb'     estimated covariance matrix for B
%
%   Example:
%       % Generate Weibull data with A depending on predictor X
%       x = 4*rand(100,1); A = 50*exp(-0.5*x); B = 2;
%       y = wblrnd(A,B);
%    
%       % Fit Cox model
%       [b,logL,H,st] = coxphfit(x,y);
%    
%       % Show Cox estimate of baseline survivor and known Weibull function
%       stairs(H(:,1),exp(-H(:,2)))
%       xx = linspace(0,100);
%       line(xx,1-wblcdf(xx,50*exp(-0.5*mean(x)),B),'color','r')
%       title(sprintf('Baseline survivor function for X=%g',mean(x)));
%
%   See also ECDF, STATSET, WBLFIT.

%   Copyright 2005-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.4 $  $Date: 2011/07/20 00:08:19 $

narginchk(2,inf);
% Check the required data arguments
if ndims(X)>2 || ~isreal(X)
    error(message('stats:coxphfit:BadX'));
end
if ~isvector(y) || ~isreal(y)
    error(message('stats:coxphfit:BadY'));
end

% Process the optional arguments
okargs =   {'baseline' 'censoring' 'frequency' 'init' 'options'};
defaults = {[]         []          []          []     []};
[baseX cens freq init options] = internal.stats.parseArgs(okargs,defaults,varargin{:});

if ~isempty(cens) && (~isvector(cens) || ~all(ismember(cens,0:1)))
    error(message('stats:coxphfit:BadCensoring'));
end
if ~isempty(freq) && (~isvector(freq) || ~isreal(freq) || any(freq<0))
    error(message('stats:coxphfit:BadFrequency'));
end
if ~isempty(baseX) && ~(isnumeric(baseX) && (isscalar(baseX) || ...
                             (isvector(baseX) && length(baseX)==size(X,2))))
    error(message('stats:coxphfit:BadBaseline'));
elseif isscalar(baseX)
    baseX = repmat(baseX,1,size(X,2));
end
% Make sure the inputs agree in size, and remove NaNs
freq(freq==0) = NaN;   % easiest way to deal with zero-frequency observations
[badin,~,y,X,cens,freq]=statremovenan(y,X,cens,freq);
if badin>0
    whichbad = {'Y' 'X' '''censoring''' '''frequency'''};
    error(message('stats:coxphfit:InputSizeMismatch', whichbad{ badin }));
end

% Sort by increasing time
[sorty,idx] = sort(y);
X = X(idx,:);
[n,p] = size(X);
if isempty(cens)
    cens = false(n,1);
else
    cens = cens(idx);
end
if isempty(freq)
    freq = ones(n,1);
else
    freq = freq(idx);
end

% Determine the observations at risk at each time
[~,atrisk] = ismember(sorty,flipud(sorty));
atrisk = length(sorty) + 1 - atrisk;     % "atrisk" used in nested function
tied = diff(sorty) == 0;
tied = [false;tied] | [tied;false];      % "tied" used in nested function

% Recenter X to make it better conditioned; does not affect coefficients
sumf = max(1, sum(freq));
if ~isempty(X)
    if isempty(baseX)
        baseX = (freq'*X) / sumf;
    end
    X = X - repmat(baseX,n,1);
end

% Try to diagnose some potential problems
if rank([ones(n,1), X]) < (p+1)
    if n>1 && any(all(diff(X,1)==0,1))
        warning(message('stats:coxphfit:DisallowedContantTerm'));
    else
        warning(message('stats:coxphfit:RankDeficient'));
    end
end

% Get starting values that will not be extreme
if isempty(init)
    stdX = sqrt((freq'*X.^2)/sumf)';
    b0 = zeros(size(stdX),class(X));
    t = (stdX ~= 0);
    b0(t) = 0.01 ./ stdX(t);
elseif isvector(init) && numel(init)==p && isreal(p) && isnumeric(p)
    b0 = init(:);
else
    error(message('stats:coxphfit:BadInit', p));
end

% Perform the fit
unboundwng = false;
if p==0
    b = zeros(0,1,class(X));
else
    % The default options include turning statsfminbx's display off.  This
    % function gives its own warning/error messages, and the caller can
    % turn display on to get the text output from statsfminbx if desired.
    options = statset(statset('coxphfit'),options);
    options = optimset(options);
    dflts = struct('DerivativeCheck','off', 'HessMult',[], ...
        'HessPattern',ones(p), 'PrecondBandWidth',Inf, ...
        'TypicalX',ones(p,1), 'MaxPCGIter',1, 'TolPCG',0.1);

    % Maximize the log-likelihood
    funfcn = {'fungradhess' 'coxphfit' @negloglike [] []};
    lastb = [];
    try
        [b, ~, ~, err, output] = ...
            statsfminbx(funfcn, b0, [], [], options, dflts, 1);
        if (err == 0)
            % Check for convergence failure
            if output.funcCount >= options.MaxFunEvals
                warning(message('stats:coxphfit:EvalLimit'));
            else
                warning(message('stats:coxphfit:IterLimit'));
            end
            unboundwng = true;
        elseif (err < 0)
            error(message('stats:coxphfit:NoSolution'));
        end

    catch ME
        if isequal(ME.identifier,'stats:statsfminbx:BadFunctionOutput') || ...
           isequal(ME.identifier,'MATLAB:eig:matrixWithNaNInf')
            warning(message('stats:coxphfit:EvalWarning'));
             b = lastb;
             unboundwng = true;
        else
            m = message('stats:coxphfit:FitError');
            throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
        end
    end
end

% Compute log likelihood at the solution, incl. baseline hazard
[logL,~,ddl,H,mlflag] = LogL(X,b,freq,cens,atrisk,tied,sorty);

% Try to diagnose a likelihood with no local maximum
if ~unboundwng && mlflag
        warning(message('stats:coxphfit:FitWarning'));
end

% Compute standard errors, etc.
if nargout>=4
    covb = cast(-inv(ddl),class(b));
    stats.covb = covb;
    varb = diag(covb);
    stats.beta = b;
    
    % In cases where we did not converge to a maximum, the 2nd derivative
    % matrix may not be negative definite.  It's likely a warning was
    % displayed and that varb has -Inf values.  We'll give infinite
    % variances for parameters affected by this.
    se = Inf(size(varb),class(b));
    se(varb>0) = sqrt(varb(varb>0));
    stats.se = se(:);  % se may be empty, we need a column vector
    stats.z = b ./ stats.se;
    stats.p = 2*normcdf(-abs(stats.z));
end

% -------------------------------------------
function [L,dl,ddl] = negloglike(b)
    % Compute negative log likelihood
    [L,dl,ddl] = LogL(X,b,freq,cens,atrisk,tied,sorty);
    L = -L;
    dl = -dl;
    ddl = -ddl;
    lastb = b;
end
end
    
function [L,dl,ddl,H,mlflag]=LogL(X,b,freq,cens,atrisk,tied,sorty)
    % Compute log likelihood L
    obsfreq = freq .* ~cens;
    Xb = X*b;
    r = exp(Xb);
    risksum = flipud(cumsum(flipud(freq.*r)));
    risksum = risksum(atrisk);
    L = obsfreq'*(Xb - log(risksum));

    if nargout>=2
        % Compute first derivative dL/db
        [n,p] = size(X);
        Xr = X .* repmat(r.*freq,1,p);
        Xrsum = flipud(cumsum(flipud(Xr)));
        Xrsum = Xrsum(atrisk,:);
        A = Xrsum ./ repmat(risksum,1,p);
        dl = obsfreq' * (X-A);
        if nargout>=5
            % Set the mlflag (monotone likelihood flag) to indicate if the
            % likelihood appears to be monotone, not at an optimum.  This
            % can happen if, at each of the sorted failure times, the
            % specified linear combination of the X's is larger than that
            % of all other observations at risk.
            if n>2
                mlflag = all(cens(1:end-1) | (diff(Xb)<0 & ~tied(1:end-1)));
            else
                mlflag = true;
            end
        end
    end
    if nargout>=3
        % Compute second derivative d2L/db2
        t1 = repmat(1:p,1,p);
        t2 = sort(t1);
        XXr = X(:,t1) .* X(:,t2) .* repmat(r.*freq,1,p^2);
        XXrsum = flipdim(cumsum(flipdim(XXr,1),1),1);
        XXrsum = XXrsum(atrisk,:,:) ./ repmat(risksum,[1,p^2]);
        ddl = reshape(-obsfreq'*XXrsum, [p,p]);
        ddl = ddl + A'*(A.*repmat(obsfreq,1,p));
    end
    if nargout>=4
        % Compute estimate of baseline cumulative hazard
        terms = obsfreq ./ risksum;
        H = [sorty, cumsum(terms)];
        H = H(obsfreq>0,:);
        if ~isempty(H)
            H = [H(1,1), 0; H];
        end
    end
    if nargout>1 && any(isnan(dl(:)))
        % If we ran in to problems in computing the derivatives and got a
        % NaN result, this will cause problems in statsfminbx.  Best to set
        % the objective function to infinity, causing statsfminbx to back
        % up.
        L = -Inf;
    end
end
