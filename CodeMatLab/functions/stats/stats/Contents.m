% Statistics Toolbox
% Version 8.1 (R2012b) 20-Jul-2012
%
% Distributions.
%  Parameter estimation.
%   betafit     - Beta parameter estimation.
%   binofit     - Binomial parameter estimation.
%   dfittool    - Distribution fitting tool.
%   evfit       - Extreme value parameter estimation.
%   expfit      - Exponential parameter estimation.
%   fitdist     - Distribution fitting.
%   gamfit      - Gamma parameter estimation.
%   gevfit      - Generalized extreme value parameter estimation.
%   gmdistribution - Gaussian mixture model estimation.
%   gpfit       - Generalized Pareto parameter estimation.
%   lognfit     - Lognormal parameter estimation.
%   mle         - Maximum likelihood estimation (MLE).
%   mlecov      - Asymptotic covariance matrix of MLE.
%   nbinfit     - Negative binomial parameter estimation.
%   normfit     - Normal parameter estimation.
%   paretotails - Empirical cdf with generalized Pareto tails.
%   poissfit    - Poisson parameter estimation.
%   raylfit     - Rayleigh parameter estimation.
%   unifit      - Uniform parameter estimation.
%   wblfit      - Weibull parameter estimation.
%
%  Probability density functions (pdf).
%   betapdf     - Beta density.
%   binopdf     - Binomial density.
%   chi2pdf     - Chi square density.
%   evpdf       - Extreme value density.
%   exppdf      - Exponential density.
%   fpdf        - F density.
%   gampdf      - Gamma density.
%   geopdf      - Geometric density.
%   gevpdf      - Generalized extreme value density.
%   gppdf       - Generalized Pareto density.
%   hygepdf     - Hypergeometric density.
%   lognpdf     - Lognormal density.
%   mnpdf       - Multinomial probability density function.
%   mvnpdf      - Multivariate normal density.
%   mvtpdf      - Multivariate t density.
%   nbinpdf     - Negative binomial density.
%   ncfpdf      - Noncentral F density.
%   nctpdf      - Noncentral t density.
%   ncx2pdf     - Noncentral Chi-square density.
%   normpdf     - Normal (Gaussian) density.
%   pdf         - Density function for a specified distribution.
%   poisspdf    - Poisson density.
%   raylpdf     - Rayleigh density.
%   tpdf        - T density.
%   unidpdf     - Discrete uniform density.
%   unifpdf     - Uniform density.
%   wblpdf      - Weibull density.
%
%  Cumulative Distribution functions (cdf).
%   betacdf     - Beta cumulative distribution function.
%   binocdf     - Binomial cumulative distribution function.
%   cdf         - Specified cumulative distribution function.
%   chi2cdf     - Chi square cumulative distribution function.
%   ecdf        - Empirical cumulative distribution function (Kaplan-Meier estimate).
%   evcdf       - Extreme value cumulative distribution function.
%   expcdf      - Exponential cumulative distribution function.
%   fcdf        - F cumulative distribution function.
%   gamcdf      - Gamma cumulative distribution function.
%   geocdf      - Geometric cumulative distribution function.
%   gevcdf      - Generalized extreme value cumulative distribution function.
%   gpcdf       - Generalized Pareto cumulative distribution function.
%   hygecdf     - Hypergeometric cumulative distribution function.
%   logncdf     - Lognormal cumulative distribution function.
%   mvncdf      - Multivariate normal cumulative distribution function.
%   mvtcdf      - Multivariate t cumulative distribution function.
%   nbincdf     - Negative binomial cumulative distribution function.
%   ncfcdf      - Noncentral F cumulative distribution function.
%   nctcdf      - Noncentral t cumulative distribution function.
%   ncx2cdf     - Noncentral Chi-square cumulative distribution function.
%   normcdf     - Normal (Gaussian) cumulative distribution function.
%   poisscdf    - Poisson cumulative distribution function.
%   raylcdf     - Rayleigh cumulative distribution function.
%   tcdf        - T cumulative distribution function.
%   unidcdf     - Discrete uniform cumulative distribution function.
%   unifcdf     - Uniform cumulative distribution function.
%   wblcdf      - Weibull cumulative distribution function.
%
%  Critical Values of Distribution functions.
%   betainv     - Beta inverse cumulative distribution function.
%   binoinv     - Binomial inverse cumulative distribution function.
%   chi2inv     - Chi square inverse cumulative distribution function.
%   evinv       - Extreme value inverse cumulative distribution function.
%   expinv      - Exponential inverse cumulative distribution function.
%   finv        - F inverse cumulative distribution function.
%   gaminv      - Gamma inverse cumulative distribution function.
%   geoinv      - Geometric inverse cumulative distribution function.
%   gevinv      - Generalized extreme value inverse cumulative distribution function.
%   gpinv       - Generalized Pareto inverse cumulative distribution function.
%   hygeinv     - Hypergeometric inverse cumulative distribution function.
%   icdf        - Specified inverse cumulative distribution function.
%   logninv     - Lognormal inverse cumulative distribution function.
%   nbininv     - Negative binomial inverse distribution function.
%   ncfinv      - Noncentral F inverse cumulative distribution function.
%   nctinv      - Noncentral t inverse cumulative distribution function.
%   ncx2inv     - Noncentral Chi-square inverse distribution function.
%   norminv     - Normal (Gaussian) inverse cumulative distribution function.
%   poissinv    - Poisson inverse cumulative distribution function.
%   raylinv     - Rayleigh inverse cumulative distribution function.
%   tinv        - T inverse cumulative distribution function.
%   unidinv     - Discrete uniform inverse cumulative distribution function.
%   unifinv     - Uniform inverse cumulative distribution function.
%   wblinv      - Weibull inverse cumulative distribution function.
%
%  Random Number Generators.
%   betarnd     - Beta random numbers.
%   binornd     - Binomial random numbers.
%   chi2rnd     - Chi square random numbers.
%   datasample  - Randomly sample from data, with or without replacement.
%   evrnd       - Extreme value random numbers.
%   exprnd      - Exponential random numbers.
%   frnd        - F random numbers.
%   gamrnd      - Gamma random numbers.
%   geornd      - Geometric random numbers.
%   gevrnd      - Generalized extreme value random numbers.
%   gprnd       - Generalized Pareto inverse random numbers.
%   hygernd     - Hypergeometric random numbers.
%   iwishrnd    - Inverse Wishart random matrix.
%   johnsrnd    - Random numbers from the Johnson system of distributions.
%   lognrnd     - Lognormal random numbers.
%   mhsample    - Metropolis-Hastings algorithm.
%   mnrnd       - Multinomial random vectors.
%   mvnrnd      - Multivariate normal random vectors.
%   mvtrnd      - Multivariate t random vectors.
%   nbinrnd     - Negative binomial random numbers.
%   ncfrnd      - Noncentral F random numbers.
%   nctrnd      - Noncentral t random numbers.
%   ncx2rnd     - Noncentral Chi-square random numbers.
%   normrnd     - Normal (Gaussian) random numbers.
%   pearsrnd    - Random numbers from the Pearson system of distributions.
%   poissrnd    - Poisson random numbers.
%   randg       - Gamma random numbers (unit scale).
%   random      - Random numbers from specified distribution.
%   randsample  - Random sample from finite population.
%   raylrnd     - Rayleigh random numbers.
%   slicesample - Slice sampling method.
%   trnd        - T random numbers.
%   unidrnd     - Discrete uniform random numbers.
%   unifrnd     - Uniform random numbers.
%   wblrnd      - Weibull random numbers.
%   wishrnd     - Wishart random matrix.
%
%  Quasi-Random Number Generators.
%   haltonset   - Halton sequence point set.       
%   qrandstream - Quasi-random stream.
%   sobolset    - Sobol sequence point set.
%
%  Statistics.
%   betastat    - Beta mean and variance.
%   binostat    - Binomial mean and variance.
%   chi2stat    - Chi square mean and variance.
%   evstat      - Extreme value mean and variance.
%   expstat     - Exponential mean and variance.
%   fstat       - F mean and variance.
%   gamstat     - Gamma mean and variance.
%   geostat     - Geometric mean and variance.
%   gevstat     - Generalized extreme value mean and variance.
%   gpstat      - Generalized Pareto inverse mean and variance.
%   hygestat    - Hypergeometric mean and variance.
%   lognstat    - Lognormal mean and variance.
%   nbinstat    - Negative binomial mean and variance.
%   ncfstat     - Noncentral F mean and variance.
%   nctstat     - Noncentral t mean and variance.
%   ncx2stat    - Noncentral Chi-square mean and variance.
%   normstat    - Normal (Gaussian) mean and variance.
%   poisstat    - Poisson mean and variance.
%   raylstat    - Rayleigh mean and variance.
%   tstat       - T mean and variance.
%   unidstat    - Discrete uniform mean and variance.
%   unifstat    - Uniform mean and variance.
%   wblstat     - Weibull mean and variance.
%
%  Likelihood functions.
%   betalike    - Negative beta log-likelihood.
%   evlike      - Negative extreme value log-likelihood.
%   explike     - Negative exponential log-likelihood.
%   gamlike     - Negative gamma log-likelihood.
%   gevlike     - Generalized extreme value log-likelihood.
%   gplike      - Generalized Pareto inverse log-likelihood.
%   lognlike    - Negative lognormal log-likelihood.
%   nbinlike    - Negative binomial log-likelihood.
%   normlike    - Negative normal likelihood.
%   wbllike     - Negative Weibull log-likelihood.
%
%  Probability distribution objects.
%   ProbDistUnivKernel - Univariate kernel smoothing distributions.
%   ProbDistUnivParam  - Univariate parametric distributions.
%
% Descriptive Statistics.
%   bootci      - Bootstrap confidence intervals.
%   bootstrp    - Bootstrap statistics.
%   corr        - Linear or rank correlation coefficient.
%   corrcoef    - Linear correlation coefficient (in MATLAB toolbox).
%   cov         - Covariance (in MATLAB toolbox).
%   crosstab    - Cross tabulation.
%   geomean     - Geometric mean.
%   grpstats    - Summary statistics by group.
%   harmmean    - Harmonic mean.
%   iqr         - Interquartile range.
%   jackknife   - Jackknife statistics.
%   kurtosis    - Kurtosis.
%   mad         - Median Absolute Deviation.
%   mean        - Sample average (in MATLAB toolbox).
%   median      - 50th percentile of a sample (in MATLAB toolbox).
%   mode        - Mode, or most frequent value in a sample (in MATLAB toolbox).
%   moment      - Moments of a sample.
%   nancov      - Covariance matrix ignoring NaNs.
%   nanmax      - Maximum ignoring NaNs.
%   nanmean     - Mean ignoring NaNs.
%   nanmedian   - Median ignoring NaNs.
%   nanmin      - Minimum ignoring NaNs.
%   nanstd      - Standard deviation ignoring NaNs.
%   nansum      - Sum ignoring NaNs.
%   nanvar      - Variance ignoring NaNs.
%   partialcorr - Linear or rank partial correlation coefficient.
%   prctile     - Percentiles.
%   quantile    - Quantiles.
%   range       - Range.
%   skewness    - Skewness.
%   std         - Standard deviation (in MATLAB toolbox).
%   tabulate    - Frequency table.
%   trimmean    - Trimmed mean.
%   var         - Variance (in MATLAB toolbox).
%
% Parametric Regression Analysis.
%  Analysis of Variance.
%   anova1      - One-way analysis of variance.
%   anova2      - Two-way analysis of variance.
%   anovan      - n-way analysis of variance.
%   aoctool     - Interactive tool for analysis of covariance.
%   friedman    - Friedman's test (nonparametric two-way anova).
%   kruskalwallis - Kruskal-Wallis test (nonparametric one-way anova).
%   manova1     - One-way multivariate analysis of variance.
%   manovacluster - Draw clusters of group means for manova1.
%   multcompare - Multiple comparisons of means and other estimates.
%
%  Linear Regression
%   GeneralizedLinearModel - Generalized linear model fitting.
%   dummyvar    - Dummy-variable coding.
%   glmfit      - Generalized linear model fitting.
%   glmval      - Evaluate fitted values for generalized linear model.
%   invpred     - Inverse prediction for simple linear regression.
%   leverage    - Regression diagnostic.
%   LinearModel - Linear least squares, robust, or stepwise regression.
%   lscov       - Ordinary, weighted, or generalized least-squares (in MATLAB toolbox).
%   lsqnonneg   - Non-negative least-squares (in MATLAB toolbox).
%   mnrfit      - Nominal or ordinal multinomial regression model fitting.
%   mnrval      - Predict values for nominal or ordinal multinomial regression.
%   multcompare - Multiple comparisons of means and other estimates.
%   mvregress   - Multivariate regression with missing data.
%   mvregresslike - Negative log-likelihood for multivariate regression.
%   polyconf    - Polynomial evaluation and confidence interval estimation.
%   polyfit     - Least-squares polynomial fitting (in MATLAB toolbox).
%   polyval     - Predicted values for polynomial functions (in MATLAB toolbox).
%   regress     - Multiple linear regression using least squares.
%   regstats    - Regression diagnostics.
%   ridge       - Ridge regression.
%   robustfit   - Robust regression model fitting.
%   stepwise    - Interactive tool for stepwise regression.
%   stepwisefit - Non-interactive stepwise regression.
%   x2fx        - Factor settings matrix (x) to design matrix (fx).
%
%  Nonlinear Regression.
%   coxphfit    - Cox proportional hazards regression.
%   nlinfit     - Nonlinear least-squares data fitting.
%   nlintool    - Interactive graphical tool for prediction in nonlinear models.
%   nlmefit     - Nonlinear mixed-effects data fitting.
%   nlmefitoutputfcn - Output function example for nlmefit and nlmefitsa.
%   nlmefitsa   - Fit nonlinear mixed effects model with stochastic EM algorithm.
%   nlpredci    - Confidence intervals for prediction.
%   nlparci     - Confidence intervals for parameters.
%   NonLinearModel - Nonlinear least-squares fitting.
%
%  Regression Plots.
%   addedvarplot - Created added-variable plot for stepwise regression.
%   nlintool    - Interactive graphical tool for prediction in nonlinear models.
%   polytool    - Interactive graph for prediction of fitted polynomials.
%   rcoplot     - Residuals case order plot.
%   robustdemo  - Interactive tool to compare robust and least squares fits.
%   rsmdemo     - Reaction simulation (DOE, RSM, nonlinear curve fitting).
%   rstool      - Multidimensional response surface visualization (RSM).
%
% Design of Experiments (DOE).
%   bbdesign    - Box-Behnken design.
%   candexch    - D-optimal design (row exchange algorithm for candidate set).
%   candgen     - Candidates set for D-optimal design generation.
%   ccdesign    - Central composite design.
%   cordexch    - D-optimal design (coordinate exchange algorithm).
%   daugment    - Augment D-optimal design.
%   dcovary     - D-optimal design with fixed covariates.
%   fracfactgen - Fractional factorial design generators.
%   ff2n        - Two-level full-factorial design.
%   fracfact    - Two-level fractional factorial design.
%   fullfact    - Mixed-level full-factorial design.
%   hadamard    - Hadamard matrices (orthogonal arrays) (in MATLAB toolbox).
%   lhsdesign   - Latin hypercube sampling design.
%   lhsnorm     - Latin hypercube multivariate normal sample.
%   rowexch     - D-optimal design (row exchange algorithm).
%
% Statistical Process Control (SPC).
%   capability  - Capability indices.
%   capaplot    - Capability plot.
%   controlchart - Shewhart control chart.
%   controlrules - Control rules (Western Electric or Nelson) for SPC data.
%   gagerr      - Gage repeatability and reproducibility (R&R) study.
%   histfit     - Histogram with superimposed normal density.
%   normspec    - Plot normal density between specification limits.
%   runstest    - Runs test for randomness.
%
% Multivariate Statistics.
%  Cluster Analysis.
%   cophenet    - Cophenetic coefficient.
%   cluster     - Construct clusters from LINKAGE output.
%   clusterdata - Construct clusters from data.
%   dendrogram  - Generate dendrogram plot.
%   gmdistribution - Gaussian mixture model estimation.
%   inconsistent - Inconsistent values of a cluster tree.
%   kmeans      - k-means clustering.
%   linkage     - Hierarchical cluster information.
%   pdist       - Pairwise distance between observations.
%   silhouette  - Silhouette plot of clustered data.
%   squareform  - Square matrix formatted distance.
%   optimalleaforder - optimal leaf ordering for hierarchical clustering.
%
%  Parametric Classification.
%   ClassificationDiscriminant - Linear discriminant analysis with regularization.
%   classify    - Linear discriminant analysis.
%   mahal       - Mahalanobis distance.
%   NaiveBayes  - Naive Bayes classification.
%
%  Dimension Reduction Techniques.
%   factoran    - Factor analysis.
%   nnmf        - Non-negative matrix factorization.
%   pca         - Principal components analysis (PCA)from raw data.
%   pcacov      - Principal components analysis (PCA) from covariance matrix.
%   pcares      - Residuals from principal components analysis (PCA).
%   rotatefactors - Rotation of FA or PCA loadings.
%
%  Copulas.
%   copulacdf   - Cumulative probability function for a copula.
%   copulafit   - Fit a parametric copula to data.
%   copulaparam - Copula parameters as a function of rank correlation.
%   copulapdf   - Probability density function for a copula.
%   copularnd   - Random vectors from a copula.
%   copulastat  - Rank correlation for a copula.
%
%  Nearest Neighbor Methods.
%   ClassificationKNN - K nearest neighbor classification.
%   createns    - Create a NeighborSearcher object for nearest neighbor search.
%   ExhaustiveSearcher - Nearest neighbor search object using exhaustive search.
%   knnsearch   - Find K nearest neighbors.
%   KDTreeSearcher - Nearest neighbor search object using kd-tree.
%   pdist2      - Pairwise distance between two sets of observations.
%   rangesearch - Find neighbors within specified radius.
%
%  Plotting.
%   andrewsplot - Andrews plot for multivariate data.
%   biplot      - Biplot of variable/factor coefficients and scores.
%   interactionplot - Interaction plot for factor effects.
%   maineffectsplot - Main effects plot for factor effects.
%   glyphplot   - Plot stars or Chernoff faces for multivariate data.
%   gplotmatrix - Matrix of scatter plots grouped by a common variable.
%   multivarichart - Multi-vari chart of factor effects.
%   parallelcoords - Parallel coordinates plot for multivariate data.
%
%  Other Multivariate Methods.
%   barttest    - Bartlett's test for dimensionality.
%   canoncorr   - Canonical correlation analysis.
%   cmdscale    - Classical multidimensional scaling.
%   mahal       - Mahalanobis distance.
%   manova1     - One-way multivariate analysis of variance.
%   mdscale     - Metric and non-metric multidimensional scaling.
%   mvregress   - Multivariate regression with missing data.
%   plsregress  - Partial least squares regression.
%   procrustes  - Procrustes analysis.
%
% Supervised Learning.
%  Decision Tree Methods.
%   ClassificationTree - Decision tree for classification.
%   classregtree - Classification and regression tree.
%   RegressionTree - Decision tree for regression.
%
%  Ensemble Methods.
%   CompactTreeBagger - Lightweight ensemble of bagged decision trees.
%   fitensemble - Fit an ensemble of learners.
%   TreeBagger  - Ensemble of bagged decision trees.
%
% Hypothesis Tests.
%   ansaribradley - Ansari-Bradley two-sample test for equal dispersions.
%   dwtest      - Durbin-Watson test for autocorrelation in linear regression.
%   linhyptest  - Linear hypothesis test on parameter estimates.
%   ranksum     - Wilcoxon rank sum test (independent samples).
%   runstest    - Runs test for randomness.
%   sampsizepwr - Sample size and power calculation for hypothesis test.
%   signrank    - Wilcoxon sign rank test (paired samples).
%   signtest    - Sign test (paired samples).
%   ttest       - One-sample and paired-sample t test.
%   ttest2      - Two-sample t test.
%   vartest     - One-sample test of variance.
%   vartest2    - Two-sample F test for equal variances.
%   vartestn    - Test for equal variances across multiple groups.
%   ztest       - Z test.
%
% Distribution Testing.
%   chi2gof     - Chi-square goodness-of-fit test.
%   jbtest      - Jarque-Bera test of normality.
%   kstest      - Kolmogorov-Smirnov test for one sample.
%   kstest2     - Kolmogorov-Smirnov test for two samples.
%   lillietest  - Lilliefors test of normality.
%
% Nonparametric Functions.
%   friedman    - Friedman's test (nonparametric two-way anova).
%   kruskalwallis - Kruskal-Wallis test (nonparametric one-way anova).
%   ksdensity   - Kernel smoothing density estimation.
%   ranksum     - Wilcoxon rank sum test (independent samples).
%   signrank    - Wilcoxon sign rank test (paired samples).
%   signtest    - Sign test (paired samples).
%
% Hidden Markov Models.
%   hmmdecode   - Calculate HMM posterior state probabilities.
%   hmmestimate - Estimate HMM parameters given state information.
%   hmmgenerate - Generate random sequence for HMM.
%   hmmtrain    - Calculate maximum likelihood estimates for HMM parameters.
%   hmmviterbi  - Calculate most probable state path for HMM sequence.
%
% Model Assessment.
%   confusionmat - Confusion matrix for classification algorithms.
%   crossval    - Loss estimate using cross-validation.
%   cvpartition - Cross-validation partition.
%   perfcurve   - ROC and other performance measures for classification algorithms.
%
% Model Selection.
%   lasso       - Lasso and elastic net linear regression.
%   lassoglm    - Lasso and elastic net generalized linear regression.
%   lassoPlot   - Lasso and elastic net plotting.
%   sequentialfs - Sequential feature selection.
%   stepwise    - Interactive tool for stepwise regression.
%   stepwisefit - Non-interactive stepwise regression.
%   relieff     - Importance of attributes (predictors) using ReliefF algorithm.
%
% Statistical Plotting.
%   andrewsplot - Andrews plot for multivariate data.
%   biplot      - Biplot of variable/factor coefficients and scores.
%   boxplot     - Boxplots of a data matrix (one per column).
%   cdfplot     - Plot of empirical cumulative distribution function.
%   ecdf        - Empirical cdf (Kaplan-Meier estimate).
%   ecdfhist    - Histogram calculated from empirical cdf.
%   fsurfht     - Interactive contour plot of a function.
%   gline       - Point, drag and click line drawing on figures.
%   glyphplot   - Plot stars or Chernoff faces for multivariate data.
%   gname       - Interactive point labeling in x-y plots.
%   gplotmatrix - Matrix of scatter plots grouped by a common variable.
%   gscatter    - Scatter plot of two variables grouped by a third.
%   hist        - Histogram (in MATLAB toolbox).
%   hist3       - Three-dimensional histogram of bivariate data.
%   ksdensity   - Kernel smoothing density estimation.
%   lsline      - Add least-square fit line to scatter plot.
%   normplot    - Normal probability plot.
%   parallelcoords - Parallel coordinates plot for multivariate data.
%   probplot    - Probability plot.
%   qqplot      - Quantile-Quantile plot.
%   refcurve    - Reference polynomial curve.
%   refline     - Reference line.
%   scatterhist - 2D scatter plot with marginal histograms.
%   surfht      - Interactive contour plot of a data grid.
%   wblplot     - Weibull probability plot.
%
% Data Objects
%   dataset     - Create datasets from workspace variables or files.
%   cell2dataset - Convert cell array to dataset array.
%   mat2dataset - Convert matrix to dataset array.
%   struct2dataset - Convert structure array to dataset array.
%   nominal     - Create arrays of nominal data.
%   ordinal     - Create arrays of ordinal data.
%
% Statistics Demos.
%   aoctool     - Interactive tool for analysis of covariance.
%   disttool    - GUI tool for exploring probability distribution functions.
%   polytool    - Interactive graph for prediction of fitted polynomials.
%   randtool    - GUI tool for generating random numbers.
%   rsmdemo     - Reaction simulation (DOE, RSM, nonlinear curve fitting).
%   robustdemo  - Interactive tool to compare robust and least squares fits.
%
% File Based I/O.
%   tblread     - Read in data in tabular format.
%   tblwrite    - Write out data in tabular format to file.
%   tdfread     - Read in text and numeric data from tab-delimited file.
%   caseread    - Read in case names.
%   casewrite   - Write out case names to file.
%   xptread     - Create a dataset array from a SAS XPORT format file.
%
% Utility Functions.
%   cholcov     - Cholesky-like decomposition for covariance matrix.
%   combnk      - Enumeration of all combinations of n objects k at a time.
%   corrcov     - Convert covariance matrix to correlation matrix.
%   groupingvariable - Help information for grouping variables.
%   grp2idx     - Convert grouping variable to indices and array of names.
%   hougen      - Prediction function for Hougen model (nonlinear example).
%   parallelstats - Help information for parallel computing options.
%   statget     - Get STATS options parameter value.
%   statset     - Set STATS options parameter value.
%   tiedrank    - Compute ranks of sample, adjusting for ties.
%   zscore      - Normalize matrix columns to mean 0, variance 1.

% Other Utility Functions.
%   cdfcalc     - Computation function for empirical cdf.
%   dfgetset    - Getting and setting dfittool parameters.
%   dfswitchyard - Invoking private functions for dfittool.
%   distchck    - Argument checking for cdf, pdf and inverse functions.
%   evnegloglike - Negative log-likelihood for the extreme value distribution.
%   iscatter    - Grouped scatter plot using integer grouping.
%   logn3fit    - 3-param lognormal fit function for cdffitdemo.
%   meansgraph  - Interactive means graph for multiple comparisons.
%   NeighborSearcher - Abstract class for nearest neighbor searching.
%   pcrsse      - SSE for principal components regression cross-validation in plspcrdemo.
%   piecewisedistribution - Parent class for paretotails.
%   ProbDist    - Abstract class for probability distributions.
%   ProbDistParametric  - Abstract class for parametric distributions.
%   ProbDistKernel - Abstract class for kernel smoothing distributions.
%   qrandset    - Parent class for haltonset and sobolset.
%   qrandstate  - State for quasi-random stream.
%   sobolstate  - State for Sobol quasi-random stream.
%   statdisptable - Display table of statistics.
%   wgtnormfit  - Weighted normal fit for customdist1demo.
%   wgtnormfit2 - Weighted normal fit based on log sigma for customdist1demo.
%
% HTML Demo Functions.
%   bayesdemo   - Bayesian analysis demo.
%   cdffitdemo  - CDF fitting demo.
%   cfitdfitdemo - Curve fitting vs. distribution fitting demo.
%   classdemo   - Classification demo.
%   clusterdemo - Cluster analysis demo.
%   cmdscaledemo - Classical multidimensional scaling demo.
%   copulademo  - Copula simulation demo.
%   creditratingdemo - Credit rating demo.
%   customdist1demo - Custom distribution fitting demo.
%   customdist2demo - Custom distribution fitting demo.
%   cvsequentialfsdemo - Feature selection and cross-validation demo.
%   factorandemo - Factor analysis demo.
%   gevdemo     - Generalized extreme value fitting demo.
%   glmdemo     - Generalized linear model demo.
%   gparetodemo - Generalized Pareto fitting demo.
%   mdscaledemo - Non-classical multidimensional scaling demo.
%   nonparametricCDFdemo - Nonparametric CDF and inverse demo.
%   mvplotdemo  - Multidimensional data plotting demo.
%   plspcrdemo  - Partial least squares regression and principal components regression demo.
%   orthoregdemo - Orthogonal regression demo.
%   samplesizedemo - Sample size calculation demo.
%   stattsdemo  - Statistics on time series demo.
%   survivaldemo - Survival data analysis demo.
%   wnlsdemo    - Weighted nonlinear least squares demo.
%   xform2lineardemo - Transforming nonlinear to linear demo.
%
% Obsolete Functions
%   betalik1    - replaced by betalike.
%   capable     - replaced by capability.
%   ewmaplot    - replaced by controlchart.
%   princomp    - replaced by pca.
%   schart      - replaced by controlchart.
%   svmsmoset   - replaced by options structure and other arguments.
%   treedisp    - replaced by classregtree/disp.
%   treefit     - replaced by classregtree.
%   treeprune   - replaced by classregtree/prune.
%   treetest    - replaced by classregtree/test.
%   treeval     - replaced by classregtree/eval.
%   weibcdf     - Weibull cdf, old parameter definitions.
%   weibfit     - Weibull fitting, old parameter definitions.
%   weibinv     - Weibull inv cdf, old parameter definitions.
%   weiblike    - Weibull likelihood, old parameter definitions.
%   weibpdf     - Weibull pdf, old parameter definitions.
%   weibplot    - Weibull prob plot, old parameter definitions.
%   weibrnd     - Weibull random numbers, old parameter definitions.
%   weibstat    - Weibull statistics, old parameter definitions.
%   xbarplot    - replaced by controlchart.

% Copyright 1993-2012 The MathWorks, Inc.
% Generated from Contents.m_template revision 1.1.8.10  $Date: 2012/05/04 00:11:47 $
