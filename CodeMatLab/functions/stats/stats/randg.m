%RANDG Gamma random numbers (unit scale).
%   Note: To generate gamma random numbers with specified shape and
%   scale parameters, you should call GAMRND.
%
%   R = RANDG returns a scalar random value chosen from a gamma
%   distribution with unit scale and shape.
%
%   R = RANDG(A) returns a matrix of random values chosen from gamma
%   distributions with unit scale.  R is the same size as A, and RANDG
%   generates each element of R using a shape parameter equal to the
%   corresponding element of A.
%
%   R = RANDG(A,M) returns an M-by-M matrix of random values chosen from
%   gamma distributions with shape parameters A.  A is either an M-by-M
%   matrix or a scalar.  If A is a scalar, RANDG uses that single shape
%   parameter value to generate all elements of R.
%
%   R = RANDG(A,M,N,...) or R = RANDG(A,[M,N,...]) returns an M-by-N-by-...
%   array of random values chosen from gamma distributions with shape
%   parameters A.  A is either an M-by-N-by-... array or a scalar.
%
%   RANDG produces pseudo-random numbers using the RAND and RANDN functions.
%   The sequence of numbers generated is determined by the settings of the
%   uniform random number generator that underlies RAND and RANDN.  Control
%   that shared random number generator using RNG.
%
%   For example, to create reproducible output from RANDG, reset the state of
%   the random number generator before calling RANDG:
%
%      rng default
%      r = randg(1,[10,1]);
%
%   See also GAMRND, RANDSTREAM, RAND, RANDN, RANDI.

%   References:
%
%      [1] Marsaglia, G. and Tsang, W.W. (2000) "A Simple Method for Generating
%          Gamma Variables", ACM Trans. Math. Soft. 26(3):363-372.

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2011/07/25 04:07:05 $
%
%   Mex file.
