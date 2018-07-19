function [distname,spec] = checkdistname(distname)
%CHECKDISTNAME Check distribution name and get information about it.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2011/05/09 01:27:51 $

if isstruct(distname)
    % It's okay to pass in the distribution structure in place of the name
    try
        spec = distname;
        distname = spec.code;
    catch me
        m = message('stats:ProbDistUnivParam:checkdistname:NotChar');
        throwAsCaller(MException(m.Identifier, '%s', getString(m)));
    end
    return
end

if ~ischar(distname) || size(distname,1)~=1
    m = message('stats:ProbDistUnivParam:checkdistname:NotChar');
    throwAsCaller(MException(m.Identifier, '%s', getString(m)));
end

% Check for abbreviations that are not just initial strings
distNames = {'extreme value' 'generalized extreme value' 'generalized pareto'...
             'negative binomial' 'discrete uniform' 'weibull'};
distAbbrevs = {'ev' 'gev' 'gp' 'nbin' 'unid' 'wbl'};
dist = lower(distname);
i = find(strcmpi(dist,distAbbrevs));
if ~isempty(i)
    dist = distNames{i};
end

% Get information about this distribution
spec = dfswitchyard('dfgetdistributions',dist);

% Should get exactly one
if numel(spec) > 1
    m = message('stats:ProbDistUnivParam:checkdistname:AmbiguousName',distname);
    throwAsCaller(MException(m.Identifier, '%s', getString(m)));
elseif isempty(spec)
    m = message('stats:ProbDistUnivParam:checkdistname:UnrecognizedName',distname);
    throwAsCaller(MException(m.Identifier, '%s', getString(m)));
end

distname = spec.code;

