%WSWOR Weighted sampling without replacement.
%   I = WSWOR(W,K) returns a vector containing a weighted random sample of
%   K unique integers selected without replacement from the set 1:N, where
%   N = LENGTH(W).  W is a vector of non-negative weights.  K must be less
%   than or equal to N.
%
%   Note that unless K == 1, the selection probabilities for the values in 1:N
%   are NOT given by W/SUM(W).  For example, if K == N, then the probability
%   that a given value in 1:N will be selected is 1.
%
%   See also RANDPERM, RANDI, DATASAMPLE.

%   Reference:
%      Wong, C.K. and M.C. Easton (1980) "An Efficient Method for Weighted
%      Sampling Without Replacement", SIAM Journal of Computing, 9(1):111-113.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2011/02/28 01:34:20 $
%
%   Mex file.
