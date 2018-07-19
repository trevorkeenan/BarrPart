function [x,cens,freq] = checkdata(spec,x,cens,freq)
%CHECKDATA Check data for fitting and return in consistent form
%    On return NaN and zero-frequency values will be squeezed out,
%    and if the distribution doesn't support censoring there will
%    be nothing in CENS or FREQ.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $  $Date: 2011/07/20 00:09:12 $

if isempty(x) || ~isvector(x) || size(x,2)~=1
    error(message('stats:ProbDistUnivParam:fit:BadData'));
elseif numel(x)<sum(~spec.prequired)
    error(message('stats:ProbDistUnivParam:fit:InsufficientData'));
end
n = numel(x);

% Check data, censoring, and frequency
bad = isnan(x);
if ~isempty(cens)
    if ~isvector(cens) || ~numel(cens)==n || (~islogical(cens) && ~all(isnan(cens) | cens==0 | cens==1))
        error(message('stats:ProbDistUnivParam:fit:BadCensoring'));
    elseif any(cens) && ~spec.censoring
        error(message('stats:ProbDistUnivParam:fit:CensoringNotAllowed', spec.code));
    else
        bad = bad | isnan(cens);
    end
end
if ~isempty(freq)
    if ~isvector(freq) || ~numel(freq)==n || ~all(isnan(freq) | (freq>=0 & freq==round(freq)))
        error(message('stats:ProbDistUnivParam:fit:BadFrequency'));
    else
        bad = bad | isnan(freq) | freq==0;
    end
end

% Remove bad values
x = x(~bad);
if ~isempty(cens)
    cens = cens(~bad);
end
if ~isempty(freq)
    freq = freq(~bad);
end
if ~spec.censoring && ~isempty(freq)
    x = expandInput(x,freq);
    freq = [];
end

% Check for inappropriate data
if ~spec.iscontinuous && any(x~=round(x))
    error(message('stats:ProbDistUnivParam:fit:NonIntegerX'));
end

if spec.closedbound(1)
    ok = min(x)>=spec.support(1);
else
    ok = min(x)>spec.support(1);
end
if ok
    if spec.closedbound(2)
        ok = max(x)<=spec.support(2);
    else
        ok = max(x)<spec.support(2);
    end
end
if ~ok
    lsym = '([';
    lsym = lsym(spec.closedbound(1)+1);
    usym = ')]';
    usym = usym(spec.closedbound(2)+1);
    error(message('stats:ProbDistUnivParam:fit:OutOfRangeX', lsym, sprintf( '%g', spec.support( 1 ) ), sprintf( '%g', spec.support( 2 ) ), usym));
end

% -----------------------------
function expanded = expandInput(input,freq)
%EXPANDDATA Expand out an input vector using element frequencies.
i = cumsum(freq);
j = zeros(1, i(end));
j(i(1:end-1)+1) = 1;
j(1) = 1;
expanded = input(cumsum(j));
