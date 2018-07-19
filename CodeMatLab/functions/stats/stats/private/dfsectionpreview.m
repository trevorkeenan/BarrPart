function imsource = dfsectionpreview(outlier, width, height)
% For use by DFITTOOL

%   $Revision: 1.1.8.3 $
%   Copyright 2003-2010 The MathWorks, Inc.

if nargin < 3
    width = 180;
    height = 180;
end

tempfigure=figure('units','pixels','position',[0 0 width height], ...
      'handlevisibility','off', ...
      'integerhandle','off', ...
      'visible','off', ...
      'paperpositionmode', 'auto', ...
      'color','w');

xlim = [0 4];
ylim = [0 4];
ax=axes('position',[.05 .05 .9 .9], ...
      'parent',tempfigure, ...
      'xtick',[],'ytick',[], ...
      'box','on', ...
      'visible','off', 'XLim',xlim,'YLim',ylim);

gr = [.9 .9 .9];
o = handle(outlier);

xlo = o.YLow;
if ~isempty(xlo)
    patch([0 1 1 0], [0 0 4 4], gr,'LineStyle','none','Parent',ax); 
end

xhi = o.YHigh;
if ~isempty(xhi)
    patch([3 4 4 3], [0 0 4 4], gr,'LineStyle','none','Parent',ax); 
end

if feature('HGUsingMATLABClasses')
    % We want the following (commented) line, but temporarily must circumvent
    % an opengl problem (cf internal problem report g658563).
    %x = print(tempfigure, '-RGBImage', '-opengl', '-r0');
    %
    % Following is the workaround until g658563 is resolved,
    % ie, embed a warning banner.
    % See also dfpreview.m, dfviewexcludepreview.m, nopreview.mat
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
    x=hardcopy(tempfigure,'-dzbuffer','-r0');
end
% give the image a black edge
x(1,:,:)=0; x(end,:,:)=0; x(:,1,:)=0; x(:,end,:)=0;
imsource=im2mis(x);

delete(tempfigure);
