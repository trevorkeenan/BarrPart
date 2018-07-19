classdef (AllowedSubclasses = {?classreg.regr.TermsRegression, ?NonLinearModel}) ParametricRegression < classreg.regr.Predictor
%ParametricRegression Fitted statistical parametric regression.
%   ParametricRegression is an abstract class representing a fitted
%   parametric regression model. You cannot create instances of this class
%   directly.  You must create a derived class by calling the fit method of
%   a derived class such as LinearModel, GeneralizedLinearModel, or
%   NonLinearModel.
%
%   See also LinearModel, GeneralizedLinearModel, NonLinearModel.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/01/12 21:36:18 $
    properties(GetAccess='public',SetAccess='protected')
%Formula - Linear or nonlinear formula describing the model.
%   The Formula property is an object that represents the the form of the
%   model. It may be a representation of the terms in a linear model, or a
%   representation of a mathematical expression in a nonlinear model.
%
%   See also fit.
        Formula = '';

%LogLikelihood - Log likelihood.
%   The LogLikelihood property is the log of the likelihood function
%   evaluated at the estimated coefficient values.
%
%   See also anova.
        LogLikelihood = NaN;

%DFE - Error degrees of freedom.
%   The DFE property is the degrees of freedom for error. DFE is the number
%   of observations minus the number of estimated coefficients.
%
%   See also anova.
        DFE = NaN; % sum(Subset) - number of estimated coefs

%SSE - Error sum of squares.
%   The SSE property is the sum of squares. SSE is the sum of the squared
%   residuals.
%
%   See also anova.
        SSE = NaN; % sum of squared errors, sum(y - yhat) over Subset

%SST - Total sum of squares.
%   The SST property is the total sum of squares. SST is the sum of the
%   squared deviations between the observed response values and their mean.
%
%   See also anova.
        SST = NaN; % total sum of squares, sum(y - mean(y)) over Subset

%SSR - Regression sum of squares.
%   The SSR property is the sum of squares explained by the regression. SSR
%   is the sum of the squared deviations between the fitted value and their
%   mean.
%
%   See also anova.
        SSR = NaN; % regression sum of squares, sum(yhat - mean(y)) over Subset

%CoefficientCovariance - Covariance matrix for coefficient estimates.
%   The CoefficientCovariance property is a matrix containing the variances
%   and covariances for the coefficient estimates.
%
%   See also Coefficients.
        CoefficientCovariance = zeros(0,0);

%CoefficientNames - Coefficient names.
%   The CoefficientNames property is a cell array of strings containing a
%   label for each coefficient. The label for the coefficient of the
%   constant term is "(Intercept)". The labels for other coefficients
%   indicate the terms that they multiply. When the term includes a
%   categorical predictor, the label indicates the level of that predictor.
%
%   See also Coefficients.
        CoefficientNames = cell(0,1);
    end
    properties(GetAccess='protected',SetAccess='protected')
        LogLikelihoodNull = NaN;
        Coefs = zeros(0,1);
    end
    properties(Dependent,GetAccess='public',SetAccess='protected')
%NumCoefficients - Number of coefficients.
%   The NumCoefficients property is the number of coefficients in the fitted
%   model. It includes coefficients that are set to zero when the model
%   terms are rank deficient.
%
%   See also Coefficients, NumEstimatedCoefficients.
        NumCoefficients

%NumEstimatedCoefficients - Number of estimated coefficients.
%   The NumEstimatedCoefficients property is the number of estimated
%   coefficients in the fitted model. It does not include coefficients that
%   are set to zero when the model terms are rank deficient. It is equal to
%   the degrees of freedom for regression.
%
%   See also Coefficients, NumCoefficients.
        NumEstimatedCoefficients

%Coefficients - Coefficient estimates and related statistics.
%   The Coefficients property is a table of coefficients. It is a dataset
%   that has one row for each coefficient and the following columns:
%      Estimate    Estimated coefficient value
%      SE          Standard error of the estimate
%      tStat       t statistic for a test that the coefficient is zero
%      pValue      p-value for the t statistic
%
%   To obtain any of these columns as a vector, index into the property
%   using dot notation. For example, in the LinearModel LM the estimated
%   coefficient vector is
%      beta = LM.Coefficients.Estimate
%
%   Use the coefTest or anova method to perform other tests on the
%   coefficients.
%
%   See also coefTest, anova.
        Coefficients

%Rsquared - Multiple correlation coefficient, or R-squared.
%   The Rsquared property provides the multiple correlation coefficient or
%   R-squared. It is a measure of the proportion of variance in the
%   response that is explained by the regression. The Rsquared property is
%   a structure with the following fields:
%      Ordinary    Ordinary (unadjusted) R-squared
%      Adjusted    R-squared adjusted for the number of coefficients
%
%   To obtain any of these values as a scalar, index into the property
%   using dot notation. For example, in the LinearModel LM the adjusted
%   R-squared value is
%      r2 = LM.Rsquared.Adjusted
%
%   See also Coefficients, anova.
        Rsquared

%ModelCriterion - AIC and other information criteria.
%   The ModelCriterion property provides AIC and other criteria used to
%   compare models. It is a structure with the following fields:
%      AIC    Akaike information criterion
%      AICc   Akaike information criterion corrected for sample size
%      BIC    Bayesian information criterion
%      CAIC   Consistent Akaike information criterion
%
%   To obtain any of these values as a scalar, index into the property
%   using dot notation. For example, in the LinearModel LM the AIC value is
%      aic = LM.ModelCriterion.AIC
%
%   See also Coefficients, anova.
        ModelCriterion % scalar struct
    end
    properties(Dependent,GetAccess='protected',SetAccess='protected')
        CoefSE
        
        % "Working" values - created and saved during fit, but subject to
        % being cleared out and recreated when needed
        y_r % reduced, i.e. y(Subset)
        w_r % reduced and normalized, i.e. Weights(Subset)/sum(Weights(Subset))
    end
    properties(Abstract,GetAccess='protected',SetAccess='protected')
        % Properties defined as abstract make it possible for different
        % subclasses that are siblings of each to define their own get/set
        % accessor methods, or even make the property dependent.  That doesn't
        % work deeper than one concrete level, but c.f. the strategy used for
        % Formula.
    end
    
    methods % get/set methods
        function p = get.NumCoefficients(model)
            p = length(model.Coefs);
        end
        function p = get.NumEstimatedCoefficients(model)
            p = model.NumObservations - model.DFE;
        end
        function table = get.Coefficients(model)
            if model.IsFitFromData
                table = tstats(model);
            else
                coefs = model.Coefs(:);
                se = sqrt(diag(model.CoefficientCovariance));
                table = dataset(coefs, se, ...
                                'VarNames',{'Value' 'SE'}, ...
                                'ObsNames',model.CoefficientNames);
            end
        end
        function rsq = get.Rsquared(model)
            rsq = get_rsquared(model);
        end
        function crit = get.ModelCriterion(model)
            crit = get_modelcriterion(model);
        end
        
        function se = get.CoefSE(model)
            se = sqrt(diag(model.CoefficientCovariance));
        end
        function y_r = get.y_r(model)
            if isempty(model.WorkingValues)
                y_r = create_y_r(model);
            else
                y_r = model.WorkingValues.y_r;
            end
        end
        function w_r = get.w_r(model)
            if isempty(model.WorkingValues)
                w_r = create_w_r(model);
            else
                w_r = model.WorkingValues.w_r;
            end
        end
    end % get/set methods
    
    methods(Access='public')
        function CI = coefCI(model,alpha) % a method, to accept alpha argument
 %coefCI Confidence intervals for coefficients.
%   CI = coefCI(M) computes 95% confidence intervals for the coefficients
%   in the regression model M. The output CI is a two-column matrix with
%   the lower confidence limits in column 1 and the upper confidence limits
%   in column 2.
%
%   CI = coefCI(M,ALPHA) computes 100*(1-ALPHA)% confidence intervals. The
%   default is ALPHA=0.05 for 95% confidence.
%
%   Example:
%       % Find 80% confidence limites for estimated coefficients
%       load carsmall
%       d = dataset(MPG,Weight);
%       d.Year = ordinal(Model_Year);
%       lm = LinearModel.fit(d,'MPG ~ Year + Weight + Weight^2')
%       confint = coefCI(lm,0.1)
%
%   See also LinearModel, GeneralizedLinearModel, NonLinearModel, linhyptest.
           if nargin < 2
                alpha = 0.05;
            end
            se = sqrt(diag(model.CoefficientCovariance));
            delta = se * tinv(1-alpha/2,model.DFE);
            CI = [(model.Coefs(:) - delta) (model.Coefs(:) + delta)];
        end
        
        % --------------------------------------------------------------------
        function [p,t,r] = coefTest(model,H,c)
%coefTest Linear hypothesis test on coefficients.
%   P = coefTest(M) computes the p-value for an F test that all
%   coefficient estimates in the regression model M are zero.
%
%   P = coefTest(M,H), with H a numeric matrix having one column for each
%   coefficient, performs an F test that H*B=0, where B represents the
%   coefficient vector.
%
%   P = coefTest(M,H,C) accepts a vector C having one value for each row
%   of H, and it performs an F test that H*B=C.
%
%   [P,F,R] = coefTest(...) also returns the F-statistic F and the rank R
%   of the matrix H. The F statistic has R degrees of freedom in the
%   numerator and M.DFE degrees of freedom in the denominator.
%
%   Example:
%      % Test to see if all coefficients in the denominator of the
%      % model expression could be the same
%      load reaction
%      myfun = 'rate~(b1*x2-x3/b5)/(1+b2*x1+b3*x2+b4*x3)';
%      nlm = NonLinearModel.fit(reactants,rate,myfun,[1 .05 .02 .1 2])
%      p = coefTest(nlm,[0 1 -1 0 0;0 0 1 -1 0])
%
%   See also LinearModel, GeneralizedLinearModel, NonLinearModel, linhyptest.

            % Test that Hb == c where b is the coefficient vector.
            nc = model.NumCoefficients;
            if nargin<2
                % Default is to test all terms = 0
                H = eye(nc);
            end
            if nargin < 3
                c = zeros(size(H,1),1);
            end
            % Outputs p = p-value, t = F statistic, r = rank of test.
            [p,t,r] = linhyptest(model.Coefs,model.CoefficientCovariance,c,H,model.DFE);
        end
    end % public
        
    methods(Hidden,Access='public') 
        function [varargout] = subsref(a,s)
            switch s(1).type
                case '()'
                    [varargout{1:nargout}] = subsref@classreg.regr.Predictor(a,s);
                case '{}'
                    [varargout{1:nargout}] = subsref@classreg.regr.Predictor(a,s);
                case '.'
                    % Look ahead so that references such as fit.Residuals.Response do not
                    % require creating all of fit.Residuals.  Let the built-in handle other
                    % properties.
                    propMethodName = s(1).subs;
                    if isequal(propMethodName,'Residuals')
                        [varargout{1:nargout}] = lookAhead(s,@a.get_residuals);
                    elseif isequal(propMethodName,'Fitted')
                        [varargout{1:nargout}] = lookAhead(s,@a.get_fitted);
                    elseif isequal(propMethodName,'Rsquared')
                        [varargout{1:nargout}] = lookAhead(s,@a.get_rsquared);
                    elseif isequal(propMethodName,'Diagnostics')
                        [varargout{1:nargout}] = lookAhead(s,@a.get_diagnostics);
                    elseif isequal(propMethodName,'ModelCriterion')
                        [varargout{1:nargout}] = lookAhead(s,@a.get_modelcriterion);
                    else
                        [varargout{1:nargout}] = subsref@classreg.regr.Predictor(a,s);
                    end
            end
        end
    end % hidden public
            

    methods(Access='protected')
        function model = ParametricRegression()
            model.PredictorTypes = 'numeric';
            model.ResponseType = 'numeric';
        end
        
        % --------------------------------------------------------------------
        function model = noFit(model,varNames,coefs,coefNames,coefCov)
            model = noFit@classreg.regr.Predictor(model,varNames);
            model.Coefs = coefs(:);
            ncoefs = length(model.Coefs);
            if isempty(coefNames)
                model.CoefficientNames = strcat({'b'},num2str((1:ncoefs)'));
            else
                if length(coefNames) ~= ncoefs
                    error(message('stats:classreg:regr:ParametricRegression:BadCoefNameSize'));
                end
                model.CoefficientNames = coefNames(:);
            end
            
            if nargin < 5
                model.CoefficientCovariance = zeros(ncoefs);
            else
                if ~isequal(size(coefCov),[ncoefs ncoefs])
                    error(message('stats:classreg:regr:ParametricRegression:BadCovarianceSize'));
                end
                [~,p] = cholcov(coefCov);
                if p > 0
                    error(message('stats:classreg:regr:ParametricRegression:BadCovarianceMatrix'));
                end
                model.CoefficientCovariance = coefCov;
            end
        end
        
        % --------------------------------------------------------------------
        function model = selectObservations(model,exclude,missing)
            if nargin < 3
                missing = [];
            end
            model = selectObservations@classreg.regr.Predictor(model,exclude,missing);
            
            % Populate the y_r and w_r fields in the WorkingValues structure
            model.WorkingValues.y_r = create_y_r(model);
            model.WorkingValues.w_r = create_w_r(model);
        end
            
        % --------------------------------------------------------------------
        function model = postFit(model)
            model = postFit@classreg.regr.Predictor(model);
            subset = model.ObservationInfo.Subset;
            resid_r = get_residuals(model,'raw');
            resid_r = resid_r(subset);
            yfit_r = predict(model);
            yfit_r = yfit_r(subset);
            w_r = model.w_r;
            sumw = sum(w_r);
            wtdymean = sum(w_r.*model.y_r) / sumw;
            model.SSE = sum(w_r .* resid_r.^2);
            model.SSR = sum(w_r .* (yfit_r - wtdymean).^2);
            model.SST = sum(w_r .* (model.y_r - wtdymean).^2);
            model.LogLikelihood = getlogLikelihood(model);
            model.LogLikelihoodNull = logLikelihoodNull(model);
        end
        
        % --------------------------------------------------------------------
        function y_r = create_y_r(model)
            subset = model.ObservationInfo.Subset;
            y = getResponse(model);
            y_r = y(subset);
        end
        function w_r = create_w_r(model)
            subset = model.ObservationInfo.Subset;
            w_r = model.ObservationInfo.Weights(subset);
        end

        function [f,p] = fTest(model)
            % F test for whole model (assumes constant term)
            ssr = model.SST - model.SSE;
            nobs = model.NumObservations;
            dfr = model.NumEstimatedCoefficients - 1;
            dfe = nobs - 1 - dfr;
            f = (ssr./dfr) / (model.SSE/dfe);
            p = fcdf(1./f,dfe,dfr); % upper tail
        end
        
        % -------------------- pass throughs to classreg.regr.modelutils -------------------
        function table = tstats(model)
            table = classreg.regr.modelutils.tstats(model.Coefs,sqrt(diag(model.CoefficientCovariance)), ...
                                      model.NumObservations,model.CoefficientNames);
        end
        function crit = get_rsquared(model,type)
            stats = struct('SSE',model.SSE, ...
                           'SST',model.SST, ...
                           'DFE',model.DFE, ...
                           'NumObservations',model.NumObservations, ...
                           'LogLikelihood',model.LogLikelihood, ...
                           'LogLikelihoodNull',model.LogLikelihoodNull);
            if nargin < 2
                crit = classreg.regr.modelutils.rsquared(stats,{'Ordinary' 'Adjusted'},true);
            else
                crit = classreg.regr.modelutils.rsquared(stats,type);
            end
        end
        function crit = get_modelcriterion(model,type)
            if nargin < 2
                crit = classreg.regr.modelutils.modelcriterion(model,'all',true);
            else
                crit = classreg.regr.modelutils.modelcriterion(model,type);
            end
        end
    end % protected
    
    methods(Abstract, Access='protected')
        model = fitter(model)% subclass-specific fitting algorithm, must fill in Coefs, CoefficientCovariance, and DFE
        L = getlogLikelihood(model)
        L0 = logLikelihoodNull(model)
    end % protected
    methods(Abstract,Hidden,Access='public') 
        v = varianceParam(model) % a consistent name for the dispersion/MSE/whatever
    end
    methods(Abstract, Static, Access='public')
        model = fit(varargin)
    end % public
end

function b = lookAhead(s,accessor)
if ~isscalar(s) && isequal(s(2).type,'.')
    subs = s(2).subs;
    if strcmp(subs,'Properties')
        b = accessor(); % a dataset array with non-variable reference
        s = s(2:end);
    else
        b = accessor(subs); % a vector
        s = s(3:end);
    end
else
    b = accessor(); % a dataset array
    s = s(2:end);
end
if ~isempty(s)
    b = subsref(b,s); 
end
end


