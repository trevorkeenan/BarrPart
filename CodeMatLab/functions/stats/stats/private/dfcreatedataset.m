function [ds, err, dsname, plotok] = dfcreatedataset(varargin)
%DFCREATEDATASET Create dfittool data set
%   DFCREATEDATASET(YEXPR,CEXPR,FEXPR,DSNAME) creates a data set named DSNAME
%   using the data obtained by evaluating the expressions entered in the gui
%   for Y, CENSORING, and FREQUENCY.
%
%   DFCREATEDATASET(YDATA,CDATA,FDATA,DSNAME,YNAME,CNAME,FNAME) creates a data
%   set named DSNAME using the y, censoring, and frequency data passed in via
%   the command line.  The "name" arguments are the names of these data
%   variables, or empty if they are not simple variable names.

%   $Revision: 1.1.8.2 $  $Date: 2011/08/17 21:29:35 $
%   Copyright 2003-2011 The MathWorks, Inc.

% Determine if we had data passed directly in from the command line
if nargin>=1 && ~ischar(varargin{1})
   fromgui = false;
else
   fromgui = true;
end

% Check for a data set name (we used to check for duplicate names
% before calling this function).

if nargin>=4
    dsname = varargin{4};
    if fromgui
        % Get the array of data sets
        dsdb = getdsdb;
        dset = down(dsdb);
        while(~isempty(dset))
            if strcmp(dset.name, dsname);
                err = 'datanamethesame';
                ds = [];
                plotok = false;
                return;
            end
            dset = right(dset);
        end
    end
else
   dsname = '';
end

n = min(length(varargin),3);
nameargs = cell(1,3);
if ~fromgui
   % Input data entered at the command line, maybe input names in args 5-7
   nameargs(:) = {''};
   valargs = cell(1,3);
   for j=1:n
      valargs{j} = varargin{j};
   end
   for j=5:nargin
      nameargs{j-4} = varargin{j};
   end
else
   % Input expressions from gui, get values by evaluating them
   nameargs = varargin(1:n);
   valargs = cell(1,0);
end

% Check input expressions or data vectors
[err,dvec,cvec,fvec,wmsg] = dfcheckselections(nameargs{:}, valargs{:});
if ~isempty(err)
   ds = [];
   plotok = false;
   return;
end
if ~isempty(wmsg)
   warndlg(wmsg,getString(message('stats:dfstrings:dlg_DistributionFittingTool')),'modal');
end

% If we're adding the first data set (maybe after old ones were
% deleted), turn off zooming
dffig = dfgetset('dffig');
if ~isequal(zoom(dffig,'getmode'),'off')
   dsdb = getdsdb;
   ds = down(dsdb);
   if isempty(ds)
       zoom(dffig,'off')
   end
end

% Make a new data set
ds = stats.dfdata(nameargs{:},dsname,dvec,cvec,fvec);
dsname = ds.name;

% Set the function type before plotting
setftype(ds,dfgetset('ftype'));
ds.conflev = dfgetset('conflev');

% Plot the histogram or other empirical curve if possible
if dfcanplotdata(ds,dffig)
   ds.plot = 1;
   ds.plotok = 1;   % this was the default when we created ds
   plotok = 'true'; % use text; boolean causes problems in java caller
else
   % Update the java dialog to show current flag state
   % ds.plot = 0;   % this was the default when we created ds
   ds.plotok = 0;
   plotok = 'false';
end

% Axis limits and the dirty flag were updated in the data set constructor,
% and the legend was updated in the updateplot method.
