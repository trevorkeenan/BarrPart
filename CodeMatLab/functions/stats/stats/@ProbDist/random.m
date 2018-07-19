function y = random(obj,varargin)
%RANDOM Random number generation. 
%   Y = RANDOM(PD) generates a random number drawn from the
%   distribution specified by the ProbDist object PD.
%
%   Y = RANDOM(PD,N) generates an N-by-N array Y of random numbers.
%
%   Y = RANDOM(PD,N,M,...) or Y=RANDOM(PD,[N,M,...]) generates an
%   N-by-M-by-... array of random numbers.
%
%   See also PROBDIST, RANDOM.

%   Copyright 2008 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2010/10/08 17:27:54 $

if isempty(obj.randfunc)
    error(message('stats:ProbDist:random:Undefined'));
else
    y = obj.randfunc(varargin{:});
end
