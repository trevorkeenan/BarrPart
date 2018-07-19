function [y,err] = mvncdf(varargin)
%MVNCDF Multivariate normal cumulative distribution function (cdf).
%   Y = MVNCDF(X) returns the cumulative probability of the multivariate
%   normal distribution with zero mean and identity covariance matrix,
%   evaluated at each row of X.  Rows of the N-by-D matrix X correspond to
%   observations or points, and columns correspond to variables or
%   coordinates.  Y is an N-by-1 vector.
%
%   Y = MVNCDF(X,MU,SIGMA) returns the cumulative probability of the
%   multivariate normal distribution with mean MU and covariance SIGMA,
%   evaluated at each row of X.  MU is a 1-by-D vector, and SIGMA is a
%   D-by-D symmetric, positive definite matrix.  MU can also be a scalar
%   value, which MVNCDF replicates to match the size of X.  If the
%   covariance matrix is diagonal, containing variances along the diagonal
%   and zero covariances off the diagonal, SIGMA may also be specified as a
%   1-by-D matrix containing just the diagonal. Pass in the empty
%   matrix for MU to use its default value when you want to only specify
%   SIGMA.
%
%   The multivariate normal cumulative probability at X is defined as the
%   probability that a random vector V, distributed as multivariate normal,
%   will fall within the semi-infinite rectangle with upper limits defined by
%   X, i.e., Pr{V(1)<=X(1), V(2)<=X(2), ... V(D)<=X(D)}.
%
%   Y = MVNCDF(XL,XU,MU,SIGMA) returns the multivariate normal cumulative
%   probability evaluated over the rectangle (hyper-rectangle for D>2) 
%   with lower and upper limits defined by XL and XU, respectively.
%
%   [Y,ERR] = MVNCDF(...) returns an estimate of the error in Y.  For
%   bivariate and trivariate distributions, MVNCDF uses adaptive quadrature on
%   a transformation of the t density, based on methods developed by Drezner
%   and Wesolowsky, and by Genz, as described in the references.  The default
%   absolute error tolerance for these cases is 1e-8.  For four or more
%   dimensions, MVNCDF uses a quasi-Monte Carlo integration algorithm based on
%   methods developed by Genz and Bretz, as described in the references.  The
%   default absolute error tolerance for these cases is 1e-4.  For
%   univariate distributions and when SIGMA is specified as a diagonal,
%   numerical integration is not used and ERR returns NaN.
%
%   [...] = MVNCDF(...,OPTIONS) specifies control parameters for the numerical
%   integration, used to compute Y when D > 1 and is not specified as a 
%   diagonal.  This argument can be created by a call to STATSET.  Choices 
%   of STATSET parameters are:
%
%         'TolFun'      - Maximum absolute error tolerance.  Default is 1e-8
%                         when D < 4, or 1e-4 when D >= 4.
%         'MaxFunEvals' - Maximum number of integrand evaluations allowed when
%                         D >= 4.  Default is 1e7.  Ignored when D < 4.
%         'Display'     - Level of display output.  Choices are 'off' (the
%                         default), 'iter', and 'final'.  Ignored when D < 4.
%
%   Example:
%
%      mu = [1 -1]; Sigma = [.9 .4; .4 .3];
%      [X1,X2] = meshgrid(linspace(-1,3,25)', linspace(-3,1,25)');
%      X = [X1(:) X2(:)];
%      p = mvncdf(X, mu, Sigma);
%      surf(X1,X2,reshape(p,25,25));
%
%   See also MVTCDF, MVNPDF, MVNRND, NORMCDF.

%   References:
%      [1] Drezner, Z. and G.O. Wesolowsky (1989) "On the Computation of the
%          Bivariate Normal Integral", J.Statist.Comput.Simul., 35:101-107.
%      [2] Drezner, Z. (1994) "Computation of the Trivariate Normal Integral",
%          Mathematics of Computation, 63:289-294.
%      [3] Genz, A. (2004) "Numerical Computation of Rectangular Bivariate
%          and Trivariate Normal and t Probabilities", Statistics and
%          Computing, 14(3):251-260.
%      [4] Genz, A. and F. Bretz (1999) "Numerical Computation of Multivariate
%          t Probabilities with Application to Power Calculation of Multiple
%          Contrasts", J.Statist.Comput.Simul., 63:361-378.
%      [5] Genz, A. and F. Bretz (2002) "Comparison of Methods for the
%          Computation of Multivariate t Probabilities", J.Comp.Graph.Stat.,
%          11(4):950-971.

%   Copyright 2005-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.4 $  $Date: 2012/05/08 20:43:05 $

% Strip off an options structure if it's there.
if isstruct(varargin{end})
    opts = statset(statset('mvncdf'), varargin{end});
    nin = nargin - 1;
else
    opts = statset('mvncdf');
    nin = nargin;
end

if nin < 1
    error(message('stats:mvncdf:TooFewInputs'));

elseif nin < 4 % MVNCDF(XU,MU,SIGMA)
    upperLimitOnly = true;
    XU = varargin{1};
    if ~ismatrix(XU)
        error(message('stats:mvncdf:InvalidDataMatrix'));
    end
    XL = -Inf(size(XU),class(XU));
    if nin > 1, mu = varargin{2}; else mu = []; end
    if nin > 2, Sigma = varargin{3}; else Sigma = []; end

else % MVNCDF(XL,XU,MU,SIGMA)
    upperLimitOnly = false;
    XL = varargin{1};
    XU = varargin{2};
    mu = varargin{3};
    Sigma = varargin{4};
    if ~ismatrix(XU) || ~isequal(size(XL),size(XU))
        error(message('stats:mvncdf:InvalidDataMatrices'));
    elseif any(any(XL > XU))
        error(message('stats:mvncdf:NonincreasingLimits'));
    end
end

% Get size of data.  Column vectors provisionally interpreted as multiple scalar data.
[n,d] = size(XU);
if d<1
    error(message('stats:mvncdf:TooFewDimensions'));
end

% Assume zero mean, data are already centered
if isempty(mu)
    XL0 = XL;
    XU0 = XU;

% Get scalar mean, and use it to center data
elseif isscalar(mu)
    XL0 = XL - mu;
    XU0 = XU - mu;

% Get vector mean, and use it to center data
elseif isvector(mu)
    [n2,d2] = size(mu);
    if d2 ~= d % has to have same number of coords as X
        error(message('stats:mvncdf:InputSizeMismatch'));
    elseif n2 == 1 % mean is a single row, rep it out to match data
        XL0 = bsxfun(@minus,XL,mu);
        XU0 = bsxfun(@minus,XU,mu);
    elseif n2 == n
        % if X and mu are column vectors and lengths match, provisionally
        % interpret this as multivariate data
        XL0 = XL - mu;
        XU0 = XU - mu;
    else % sizes don't match
        error(message('stats:mvncdf:InputSizeMismatch'));
    end
else
    error(message('stats:mvncdf:InvalidMu'));
end

% Assume identity covariance, data are already standardized
if isempty(Sigma)
    % Special case: if Sigma isn't supplied, then interpret X
    % and MU as row vectors if they were both column vectors
    if d == 1
        XL0 = XL0';
        XU0 = XU0';
        [n,d] = size(XU0);
    end
    sigmaIsDiag = true;
    Sigma = ones(1,d);
else
    sz = size(Sigma);
    if sz(1)==1 && sz(2)>1
        % Just the diagonal of Sigma has been passed in.
        sz(1) = sz(2);
        sigmaIsDiag = true;
    else
        sigmaIsDiag = false;
    end

    % Special case: if Sigma is supplied, then use it to try to interpret
    % X and MU as row vectors if they were both column vectors.
    if (d == 1)
        if sz(1) == n
            XL0 = XL0';
            XU0 = XU0';
            [n,d] = size(XU0);
        elseif ~isscalar(mu)
            error(message('stats:mvncdf:InputSizeMismatch'));
        end
    end

    % Make sure Sigma is a valid covariance matrix
    if sz(1) ~= sz(2)
        error(message('stats:mvncdf:BadCovariance'));
    elseif ~isequal(sz, [d d])
        error(message('stats:mvncdf:CovSizeMismatch'));
    else
        if ~sigmaIsDiag
            [~,err] = cholcov(Sigma,0);
            if err ~= 0
                error(message('stats:mvncdf:BadMatrixSigma'));
            end
        else
            if any(Sigma<=0)
                error(message('stats:mvncdf:BadDiagSigma'));
            end
        end
    end
end

% Standardize Sigma and X to correlation if necessary
if sigmaIsDiag
    t = sqrt(Sigma);
    % Rho is not used by special case formula used for diagonal sigma.
    %Rho = eye(d); 
    XL0 = bsxfun(@rdivide,XL0,t);
    XU0 = bsxfun(@rdivide,XU0,t);
else
    s = sqrt(diag(Sigma));
    Rho = Sigma ./ (s*s');
    XL0 = bsxfun(@rdivide,XL0,s');
    XU0 = bsxfun(@rdivide,XU0,s');
end

% Call the appropriate routine to compute the cdf from standardized values.
if d == 1 
    y = normcdf(XU0,0,1) - normcdf(XL0,0,1);
    if nargout > 1
        err = NaN(size(y),class(y));
    end
    
elseif sigmaIsDiag
    y = prod(normcdf(XU0,0,1)-normcdf(XL0,0,1),2);
    if nargout > 1
        err = NaN(size(y),class(y));
    end
    
elseif d <= 3
    tol = opts.TolFun; if isempty(tol), tol = 1e-8; end
    if d == 2, rho = Rho(2); else rho = Rho([2 3 6]); end
    if upperLimitOnly
        if d == 2
            y = bvncdf(XU0, rho, tol);
        else
            y = tvncdf(XU0, rho, tol);
        end
    else % lower and upper limits
        % Compute the probability over the rectangle as sums and differences
        % of integrals over semi-infinite half-rectangles.  For degenerate
        % rectangles, force an exact zero by making each piece exactly zero.
        equalLims = (XL0==XU0);
        XL0(equalLims) = -Inf;
        XU0(equalLims) = -Inf;
        y = zeros(n,1,superiorfloat(XL0,XU0,Rho));
        for i = 0:d
            k = nchoosek(1:d,i);
            for j = 1:size(k,1)
                X = XU0; X(:,k(j,:)) = XL0(:,k(j,:));
                if d == 2
                    y = y + (-1)^i * bvncdf(X, rho, tol/4);
                else
                    y = y + (-1)^i * tvncdf(X, rho, tol/8);
                end
            end
        end
    end
    if nargout > 1
        err = repmat(cast(tol,class(y)),size(y));
    end

elseif d <= 25
    tol = opts.TolFun; if isempty(tol), tol = 1e-4; end
    maxfunevals = opts.MaxFunEvals;
    verbose = find(strcmp(opts.Display,{'off' 'final' 'iter'})) - 1;
    y = zeros(n,1,superiorfloat(XL0,XU0,Rho));
    err = zeros(n,1,class(y));
    for i = 1:n
        [y(i),err(i)] = mvtcdfqmc(XL0(i,:),XU0(i,:),Rho,Inf,tol,maxfunevals,verbose);
    end

else
    error(message('stats:mvncdf:DimensionTooLarge'));
end

% repair roundoff or subtraction cancellation problems; min and max would drop NaNs
y(y<0) = 0;
y(y>1) = 1; 
end


%----------------------------------------------------------------------
function p = bvncdf(b,rho,~)
% CDF for the bivariate normal.
%
% Implements Section 2.4 of Genz (2004).

n = size(b,1);
if rho == 0
    p = cast(prod(Phi(b),2), superiorfloat(b,rho));
else
    if abs(rho) < 0.3      % 6 point Gauss Legendre abscissas and weights
        w = [0.4679139345726904  0.3607615730481384  0.1713244923791705];
        y = [0.2386191860831970  0.6612093864662647  0.9324695142031522];
    elseif abs(rho) < 0.75 % 12 point Gauss Legendre abscissas and weights
        w = [0.2491470458134029  0.2334925365383547  0.2031674267230659 ...
             0.1600783285433464  0.1069393259953183  0.04717533638651177];
        y = [0.1252334085114692  0.3678314989981802  0.5873179542866171 ...
             0.7699026741943050  0.9041172563704750  0.9815606342467191];
    else                 % 20 point Gauss Legendre abscissas and weights
        w = [0.1527533871307259  0.1491729864726037  0.1420961093183821  0.1316886384491766  0.1181945319615184 ...
             0.1019301198172404  0.08327674157670475 0.06267204833410906 0.04060142980038694 0.01761400713915212];
        y = [0.07652652113349733 0.2277858511416451  0.3737060887154196  0.5108670019508271  0.6360536807265150 ...
             0.7463319064601508  0.8391169718222188  0.9122344282513259  0.9639719272779138  0.9931285991850949];
    end
    
    if abs(rho) < .925
        p1 = prod(Phi(b),2);
        asinrho = asin(rho);
        w = fliplr([w w]);
        theta = asinrho .* (fliplr([y -y]) + 1) / 2;
        sintheta = sin(theta);
        cossqtheta = cos(theta).^2; % always positive
        h = -b(:,1); k = -b(:,2);
        hk = h .* k;
        ssq = (h.^2 + k.^2) / 2;
        p2 = zeros(size(p1),class(rho));
        for i = 1:n
            if isfinite(hk(i))
                f = exp( -(ssq(i) - hk(i)*sintheta) ./ cossqtheta );
                p2(i) = asinrho * sum(w.*f) / 2;
            else
                % This piece is zero if either limit is +/- infinity.  If
                % either is NaN, p1 will already be NaN.
            end
        end
        p = p1 + p2/(2*pi);

    else % abs(rho) >= .925
        if rho > 0
            p1 = Phi(min(b,[],2));
            p1(any(isnan(b),2)) = NaN;
        else
            p1 = Phi(b(:,1)) - Phi(-b(:,2));
            p1(p1<0) = 0; % max would drop NaNs
        end
        
        s = sign(rho);
        if abs(rho) < 1
            h = -b(:,1); k = -b(:,2);
            shk = s.*h.*k;
            asq = 1 - rho.^2;  a = sqrt(asq);
            b = abs(h - s.*k); bsq = b.^2;
            c = (4 - shk)/8;
            d = c .* (12 - shk)/16;

            t1 = a.*(1 + d.*asq.^2/5 + (c/3 - d.*bsq/15).*(asq - bsq));
            t2 = b.*(1 - c.*bsq/3 + d.*bsq.^2/15);
            p2 = exp(-shk/2) ...
                .* (t1.*exp(-bsq./(2*asq)) - t2.*sqrt(2*pi).*Phi(-b./a));

            w = [w w];
            x = a .* ([y -y] + 1) / 2; x2 = x.^2; x4 = x2.^2;
            sqrt1mx2 = sqrt(1 - x2);
            t = (sqrt1mx2 - 1) ./ (2*(sqrt1mx2 + 1));
            for i = 1:n
                if isfinite(shk(i))
                    f = exp(-(bsq(i)./x2 + shk(i))/2) ...
                        .* (exp(shk(i).*t)./sqrt1mx2 - (1 + c(i).*x2 + d(i).*x4));
                    p2(i) = p2(i) + a .* sum(w.*f) / 2;
                else
                    % This piece is zero if either limit is +/- infinity.  If
                    % either is NaN, p1 will already be NaN.
                    p2(i) = 0;
                end
            end
        else % abs(rho) == 1
            p2 = zeros(class(rho));
        end
        p = p1 - s.*p2/(2*pi);
    end
end

end % bvncdf


%----------------------------------------------------------------------
function p = tvncdf(b,rho,tol)
% CDF for the trivariate normal
%
% Implements equation (14) in Section 3.2 of Genz (2004), integrating each
% term in (14) separately in terms of theta between 0 and asin(rho_j1), using
% adaptive quadrature.

n = size(b,1);

% Find a permutation that makes rho_32 == max(rho)
[dum,imax] = max(abs(rho)); %#ok<ASGLU>
if imax == 1 % swap 1 and 3
    rho_21 = rho(3); rho_31 = rho(2); rho_32 = rho(1);
    b = b(:,[3 2 1]);
elseif imax == 2 % swap 1 and 2
    rho_21 = rho(1); rho_31 = rho(3); rho_32 = rho(2);
    b = b(:,[2 1 3]);
else % imax == 3
    rho_21 = rho(1); rho_31 = rho(2); rho_32 = rho(3);
    % b already in correct order
end

p1 = Phi(b(:,1)).*bvncdf(b(:,2:3),rho_32,tol/3);

if abs(rho_21) > 0
    loLimit = 0;
    hiLimit = asin(rho_21);
    rho_j1 = rho_21;
    rho_k1 = rho_31;
    p2 = zeros(size(p1),class(p1));
    for i = 1:n
        b1 = b(i,1); bj = b(i,2); bk = b(i,3);
        if isfinite(b1) && isfinite(bj) && ~isnan(bk)
            p2(i) = quadgk(@tvnIntegrand,loLimit,hiLimit,'AbsTol',tol/3,'RelTol',0);
        else
            % This piece is zero if either limit is +/- infinity.  If
            % either is NaN, p1 will already be NaN.
        end
    end
else
    p2 = zeros(class(p1));
end

if abs(rho_31) > 0
    loLimit = 0;
    hiLimit = asin(rho_31);
    rho_j1 = rho_31;
    rho_k1 = rho_21;
    p3 = zeros(size(p1),class(p1));
    for i = 1:n
        b1 = b(i,1); bj = b(i,3); bk = b(i,2);
        if isfinite(b1) && isfinite(bj) && ~isnan(bk)
            p3(i) = quadgk(@tvnIntegrand,loLimit,hiLimit,'AbsTol',tol/3,'RelTol',0);
        else
            % This piece is zero if either limit is +/- infinity.  If
            % either is NaN, p1 will already be NaN.
        end
    end
else
    p3 = zeros(class(p1));
end

p = cast(p1 + (p2 + p3)./(2.*pi), superiorfloat(b,rho));

    function integrand = tvnIntegrand(theta)
        % Integrand is exp( -(b1.^2 + bj.^2 - 2*b1*bj*sin(theta))/(2*cos(theta).^2) )
        sintheta = sin(theta);
        cossqtheta = cos(theta).^2; % always positive
        expon = ((b1*sintheta - bj).^2 ./ cossqtheta + b1.^2)/2;

        sinphi = sintheta .* rho_k1 ./ rho_j1;
        numeru = bk.*cossqtheta - b1.*(sinphi - rho_32.*sintheta) ...
                                - bj.*(rho_32 - sintheta.*sinphi);
        denomu = sqrt(cossqtheta.*(cossqtheta -sinphi.*sinphi ...
                                   - rho_32.*(rho_32 - 2.*sintheta.*sinphi)));
        integrand = exp(-expon) .* Phi(numeru./denomu);
    end
end % tvncdf


%----------------------------------------------------------------------
function p = Phi(z)
% CDF for the normal distribution.
p = 0.5 * erfc(-z / sqrt(2));
end
