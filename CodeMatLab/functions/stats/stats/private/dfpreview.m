function [err, imsource] = dfpreview(dexpr, cexpr, fexpr, width, height, ds, binInfo)
% For use by DFITTOOL

%   $Revision: 1.1.8.4 $  $Date: 2010/11/08 02:37:43 $
%   Copyright 2001-2010 The MathWorks, Inc.

if nargin<6
    [err, data, censoring, frequency] = dfcheckselections(dexpr, cexpr, fexpr);
    if ~isequal(err, '')
        imsource = [];
        return;
    end
else
    err = '';
    data = ds.y;
    censoring = ds.censored;
    frequency=ds.frequency;
end

if nargin < 5
    width = 200;
    height = 200;
end

tempfigure=figure('units','pixels','position',[0 0 width height], ...
    'handlevisibility','callback', ...
    'integerhandle','off', ...
    'visible','off', ...
    'paperpositionmode', 'auto', ...
    'color','w');
tempaxes=axes('position',[.05 .05 .9 .9], ...
   'parent',tempfigure, ...
   'box','on', ...
   'visible','off');

% If data has a complex part, it will spit a warning to the command line, so
% turn off warnings before plotting.
warnstate=warning('off', 'all');

if nargin < 6
    binInfo = dfgetset('binDlgInfo');
elseif nargin < 7 
    binInfo = ds.binDlgInfo;
else
    % binInfo passed in
end

% If we're working on expressions rather than data in an existing data set,
% we may need to remove NaNs
[ignore1,ignore2,data,censoring,frequency] = statremovenan(data,censoring,frequency);

% Compute the bin centers using the ecdf
% to allow a quartile computation even when there is censoring.
[fstep, xstep] = ecdf(data, 'censoring', censoring, 'frequency', frequency);
[dum,binEdges] = internal.stats.histbins(data,censoring,frequency,binInfo,fstep,xstep);

set(0,'CurrentFigure', tempfigure);
set(tempfigure,'CurrentAxes', tempaxes);

% Plot a histogram from ecdf using the computed number of bins
ecdfhist(tempaxes, fstep, xstep, 'edges', binEdges);
set(tempaxes, 'xtick',[],'ytick',[]);
axis(tempaxes,'tight');
allchildren = get(tempaxes, 'children');
patchchildren = findobj(allchildren,'flat','Type','patch');
set(patchchildren, 'facecolor', [.9 .9 .9]);
warning(warnstate);

if feature('HGUsingMATLABClasses')
    % We want the following (commented) line, but temporarily must circumvent
    % an opengl problem (cf internal problem report g658563).
    %x = print(tempfigure, '-RGBImage', '-opengl', '-r0');
    %
    % Following is the workaround until g658563 is resolved,
    % ie, embed a warning banner.
    % See also dfviewexcludepreview.m, dfsectionpreview.m, nopreview.mat
    x = ones(height,width,3,'uint8') * uint8(255);   
    banner = load('nopreview'); 
    sz = size(banner.img); 
    clipht = min(height,sz(1));
    clipwd = min(width,sz(2));
    hoffset = floor((height-sz(1))/2); 
    woffset = floor((width-sz(2))/2); 
    bhoffset = 0;
    bwoffset = 0;
    if hoffset < 0
        bhoffset = -hoffset;
        hoffset = 0;
    end
    if woffset < 0
        bwoffset = -woffset;
        woffset = 0;
    end
    x(hoffset+1:hoffset+clipht, woffset+1:woffset+clipwd, :) = ...
        banner.img(bhoffset+1:bhoffset+clipht, bwoffset+1:bwoffset+clipwd, :); 
else
    x = hardcopy(tempfigure,'-dzbuffer','-r0');
end
% give the image a black edge
x(1,:,:)=0; x(end,:,:)=0; x(:,1,:)=0; x(:,end,:)=0;
imsource=im2mis(x);

delete(tempfigure);
