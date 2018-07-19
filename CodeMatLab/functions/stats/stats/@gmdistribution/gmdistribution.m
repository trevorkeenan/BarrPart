classdef gmdistribution
%GMDISTRIBUTION Gaussian mixture distribution class.
%   An object of the GMDISTRIBUTION class defines a Gaussian mixture
%   distribution, which is a multivariate distribution that consists of a
%   mixture of one or more multivariate Gaussian distribution components.  The
%   number of components for a given GMDISTRIBUTION object is fixed.  Each
%   multivariate Gaussian component is defined by its mean and covariance, and
%   the mixture is defined by a vector of mixing proportions.
%
%   To create a Gaussian mixture distribution by specifying the distribution
%   parameters, use the GMDISTRIBUTION constructor.  To fit a Gaussian mixture
%   distribution model to data, use GMDISTRIBUTION.FIT.
%
%   A Gaussian mixture distribution with K components, in D dimensions, has
%   the following properties:
%
%      NDimensions  The number of dimensions for each of the multivariate
%                   Gaussian components in the mixture distribution, D.
%      DistName     'gaussian mixture distribution', the name of the
%                   distribution.
%      NComponents  The number of mixture components, K.
%      PComponents  A 1-by-K vector containing the mixing proportion of
%                   each component.
%      CovType      'diagonal' if the component covariance matrices are
%                   restricted to be diagonal; 'full' otherwise. 
%      SharedCov    True if all the component covariance matrices are
%                   restricted to be the same (pooled cov); false otherwise.
%      mu           A K-by-D matrix of component means.
%      Sigma        An array or a matrix containing the component covariance
%                   matrices.  Sigma is one of the following
%                      * A D-by-D-by-K array if there are no restrictions on
%                        the form of covariance.  In this case, Sigma(:,:,J)
%                        is the covariance matrix of component J.
%                      * A 1-by-D-by-K array if the covariance matrices are
%                        restricted to be diagonal, but not restricted to be
%                        same across components.  In this case Sigma(:,:,J)
%                        contains the diagonal elements of the covariance
%                        matrix of component J.
%                      * A D-by-D matrix if the covariance matrices are
%                        restricted to be the same across components, but not
%                        restricted to be diagonal.  In this case, Sigma is
%                        the common covariance matrix.
%                      * A 1-by-D vector if the covariance matrices are
%                        restricted to be diagonal and to be the same across
%                        components.  In this case, Sigma contains the
%                        diagonal elements of the common covariance matrix.
%
%   A Gaussian mixture distribution object created by fitting to data using
%   GMDISTRIBUTION.FIT also has the following properties:
%
%      NlogL        The negative of the log-likelihood of the fit.
%      AIC          The Akaike information criterion for the fit, defined as
%                   2*NlogL + 2*(the number of estimated parameters).
%      BIC          The Bayes information criterion for the fit, defined as
%                   2*NlogL + (the number of estimated parameters * log(N)).             
%      Converged    True if the fitting algorithm converged; false if the
%                   algorithm did not converge.
%      Iters        The number of iterations taken by the fitting algorithm.
%      RegV         The value supplied for the 'Regularize' input parameter
%                   to the FIT method.
%
%   See also GMDISTRIBUTION/GMDISTRIBUTION, GMDISTRIBUTION/FIT,
%            GMDISTRIBUTION/CLUSTER, GMDISTRIBUTION/PDF, GMDISTRIBUTION/CDF,
%            GMDISTRIBUTION/RANDOM, GMDISTRIBUTION/POSTERIOR,
%            GMDISTRIBUTION/MAHAL. 

%   Reference:   McLachlan, G., and D. Peel, Finite Mixture Models, John
%                Wiley & Sons, New York, 2000.

%   Copyright 2007-2010 The MathWorks, Inc.
%   $Revision: 1.1.8.5 $  $Date: 2011/05/09 01:28:07 $

    properties(GetAccess='public', SetAccess='protected')
        NDimensions = 0;  % dimension of multivariate distribution
        DistName = 'gaussian mixture distribution';
        NComponents = 0;  % number of mixture components
        PComponents = zeros(0,1); % NComponents-by-1 vector of proportions
        mu = [];          % NComponents-by-NDimensions array for means
        Sigma = [];       % Covariance
        NlogL=[];         % Negative log-likelihood
        AIC = [];         % Akaike information criterion
        BIC = [];         % Bayes information criterion
        Converged = [];   % Has the EM converged
        Iters = [];       % The number of iterations
        SharedCov = [];
        CovType = [];
        RegV = 0;
    end

    methods
        function obj = gmdistribution(mu, Sigma, p)
        %GMDISTRIBUTION Create a Gaussian mixture model.
        %   GM = GMDISTRIBUTION(MU,SIGMA,P) creates a distribution consisting
        %   of a mixture of multivariate Gaussian components, given values for
        %   the components' distribution parameters.  To create a Gaussian
        %   mixture distribution by fitting to data, use GMDISTTRIBUTION.FIT.
        %
        %   The number of components and the dimension of the distribution are
        %   implicitly defined by the sizes of the inputs MU, SIGMA, and P.
        %
        %   MU is K-by-D matrix specifying the mean of each component, where K
        %   is the number of components, and D is the number of dimensions.  MU(J,:)
        %   is the mean of component J.
        %
        %   SIGMA specifies the covariance matrix of each component.  SIGMA is one
        %   of the following:
        %
        %      * A D-by-D-by-K array if there are no restrictions on the form of the
        %        covariance matrices.  In this case, SIGMA(:,:,J) is the covariance
        %        matrix of component J.
        %      * A 1-by-D-by-K array if the covariance matrices are restricted to be
        %        diagonal, but not restricted to be same across components.  In this
        %        case, SIGMA(:,:,J) contains the diagonal elements of the covariance
        %        matrix of component J.
        %      * A D-by-D matrix if the covariance matrices are restricted to be the
        %        same across components, but not restricted to be diagonal.  In this
        %        case, SIGMA is the common covariance matrix.
        %      * A 1-by-D vector if the covariance matrices are restricted to be
        %        diagonal and the same across components.  In this case, SIGMA contains
        %        the diagonal elements of the common covariance matrix.
        %
        %   P is 1-by-K vector specifying the mixing proportions of each component.  If
        %   P does not sum to 1, GMDISTRIBUTION normalizes it.  The default is equal
        %   proportions if P is not given.
        %
        %   The inputs MU, SIGMA, and P are stored in the mu, Sigma, and PComponents
        %   properties, respectively, of GM.
        %
        %   Example:  Create a 2-component Gaussian mixture model.
        %
        %            mu = [1 2;-3 -5];
        %            Sigma = cat(3,[2 0; 0 .5],[1 0; 0 1]);
        %            mixp = ones(1,2)/2;
        %            gm = gmdistribution(mu,Sigma,mixp);
        %
        %   See also GMDISTRIBUTION, GMDISTRIBUTION/FIT.
        
            if nargin==0
                return;
            end

            if nargin < 2
                error(message('stats:gmdistribution:TooFewInputs'));
            end
            if ndims(mu) ~= 2 || ~isnumeric(mu)
                error(message('stats:gmdistribution:BadMu'));
            end

            [k,d] = size(mu);
            if nargin < 3 || isempty(p)
                p = ones(1,k);
            elseif ~isvector(p) || length(p) ~= k
                error(message('stats:gmdistribution:MisshapedMuP'));
                     
            elseif any(p <= 0)
                error(message('stats:gmdistribution:InvalidP'));
                      
            elseif size(p,1) ~= 1
                p = p'; % make it a row vector
            end

            p = p/sum(p);

            [d1,d2,k2] = size(Sigma);
            if  k2 == 1
                if d1 == 1 %diagonal covariance
                    if d2 ~= d
                        error(message('stats:gmdistribution:MisshapedSharedDiagCov'));
                    elseif any(Sigma<0)
                        error(message('stats:gmdistribution:BadDiagCov'));
                    end
                    obj.CovType = 'diagonal';

                else %full covariance
                    if ~isequal(size(Sigma),[d d])
                        error(message('stats:gmdistribution:MisshapedSharedCov'));
                    end
                    [~,err] = cholcov(Sigma);
                    if err ~= 0
                        error(message('stats:gmdistribution:BadSharedCov'));
                    end
                    obj.CovType = 'full';
                end

                obj.SharedCov = true;
                
            else %different covariance
                if k2 ~= k
                    error(message('stats:gmdistribution:MisshapedMuCov'));
                end
                if d1 == 1 %diagonal covariance
                    if d2 ~= d
                        error(message('stats:gmdistribution:MisshapedDiagCov'));
                     end
                    for j = 1:k
                        %check whether the covariance matrix is positive definite
                        if any(Sigma(:,:,j)<0)
                            error(message('stats:gmdistribution:BadDiagCov'));
                        end
                    end
                    obj.CovType = 'diagonal';
                else
                    if d1 ~= d || d2 ~= d
                        error(message('stats:gmdistribution:MisshapedCov'));
                    end
                    for j = 1:k
                        % Make sure Sigma is a valid covariance matrix
                        % check for positive definite
                        [~,err] = cholcov(Sigma(:,:,j));
                        if err ~= 0
                            error(message('stats:gmdistribution:BadCov'));
                        end
                    end
                    obj.CovType = 'full';
                end
                obj.SharedCov = false;
            end

            obj.NDimensions = d;
            obj.NComponents = k;
            obj.PComponents = p;
            obj.mu = mu;
            obj.Sigma = Sigma;
        end % constructor
    end

    methods(Static = true)
        obj = fit(X,k,varargin);
    end

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
            throwUndefinedError();
        end
    end
   
end % classdef

function throwNoCatError()
error(message('stats:gmdistribution:NoCatAllowed',upper(mfilename))); 
end

function throwUndefinedError()
st = dbstack;
name = regexp(st(2).name,'\.','split');
error(message('stats:gmdistribution:UndefinedFunction', name{ 2 }, mfilename));
end

