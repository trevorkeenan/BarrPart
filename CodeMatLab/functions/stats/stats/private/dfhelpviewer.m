function dfhelpviewer(topic, errorname)
% DFHELPVIEWER  is a helper file for the Distribution Fitting Toolbox 
% DFHELPVIEWER Displays help for Distribution Fitting TOPIC. If the map file 
% cannot be found, an error is displayed using ERRORNAME

%   Copyright 2003-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $

mapfilename = [docroot '/toolbox/stats/stats.map'];
try
    helpview(mapfilename, topic);
catch
    errordlg(getString(message('stats:dfstrings:sprintf_UnableToDisplayHelp', errorname)));
end
