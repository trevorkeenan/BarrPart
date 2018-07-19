function x = betainv(p,a,b)
%BETAINV Inverse of the beta cumulative distribution function (cdf).
%   X = BETAINV(P,A,B) returns the inverse of the beta cdf with 
%   parameters A and B at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   BETAINV uses Newton's method to converge to the solution.
%
%   See also BETACDF, BETAFIT, BETALIKE, BETAPDF, BETARND, BETASTAT, ICDF.

%   Reference:
%      [1]     M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964.

%   Copyright 1993-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2011/05/09 01:24:09 $

if nargin < 3
    error(message('stats:betainv:TooFewInputs'));
end

[errorcode, p, a, b] = distchck(3,p,a,b);
if errorcode > 0
    error(message('stats:betainv:InputSizeMismatch'));
end

% Weed out any out of range parameters or probabilities.
okAB = (0 < a & a < Inf) & (0 < b & b < Inf);
k = (okAB & (0 <= p & p <= 1));
allOK = all(k(:));

% Fill in NaNs for out of range cases.
if ~allOK
    if isa(p,'single') || isa(a,'single') || isa(b,'single')
       x = NaN(size(k),'single');
    else
       x = NaN(size(k));
    end

    % Remove the out of range cases.  If there's nothing remaining, return.
    if any(k(:))
        if numel(p) > 1, p = p(k); end
        if numel(a) > 1, a = a(k); end
        if numel(b) > 1, b = b(k); end
    else
        return;
    end
end

% Call BETAINCINV to find a root of BETAINC(Q,A,B) = P
q = betaincinv(p,a,b);

badcdf = ((abs(betainc(q,a,b) - p)./p) > sqrt(eps(class(q))));
if any(badcdf(:))   % cdf is too far off
    didnt = find(badcdf, 1, 'first');
    if numel(a) == 1, abad = a; else abad = a(didnt); end
    if numel(b) == 1, bbad = b; else bbad = b(didnt); end
    if numel(p) == 1, pbad = p; else pbad = p(didnt); end
    warning(message('stats:betainv:NoConvergence', sprintf( '%g', abad ), sprintf( '%g', bbad ), sprintf( '%g', pbad )));
end

% Broadcast the values to the correct place if need be.
if allOK
    x = q;
else
    x(k) = q;
end
