function dfsave2ws(fitobj)
%DFSAVE2WS Utility to save probability distribution object from dfittool

%   $Revision: 1.1.6.4 $  $Date: 2011/08/17 21:29:46 $
%   Copyright 2003-2011 The MathWorks, Inc.

if ischar(fitobj)
    % Fit name passed in, so get fit object
    fitdb = getfitdb;
    fitobj = find(fitdb, 'name', fitobj);
end

if ~isa(fitobj,'stats.dffit')
    error(message('stats:dfsave2w:BadFit'));
end

% Bring up dialog to get variable name for this fit
export2wsdlg({getString(message('stats:dfstrings:dlg_SaveFittedDistributionAs'))},...
             {'pd'},{fitobj.probdist},...
             getString(message('stats:dfstrings:dlg_SaveFit2WS')));

end
