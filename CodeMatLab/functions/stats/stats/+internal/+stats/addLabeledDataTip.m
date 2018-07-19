function addLabeledDataTip(labels,hLabel,hUnlabel,linelabels)
%addLabeledDataTip Define data tip for handle based on labels.
%    addLabeledDataTip(LABELS,HLABEL,HUNLABEL) accepts a cell array of
%    strings LABELS to be used in the data tip for the handles in the array
%    HLABEL. Handles listed in neither HLABEL nor HUNLABEL have the
%    default data tips. Handles in the array HUNLABEL are specified as not
%    having data tips; these may be lines such as grid or reference lines
%    that do not represent data values.
%
%    Alternatively, LABELS can be:
%        empty    to use 'Observation N' as the label for point N
%        fun      to use the result of fun(N)
%
%    addLabeledDataTip(LABELS,HLABEL,HUNLABEL,LINELABELS) also accepts an
%    optional cell array of strings of the same length as HLABEL. These
%    provide an additional label to appear at the top of the data tip for
%    the corresponding line.

%   $Revision: 1.1.6.3 $  $Date: 2012/03/01 02:30:03 $
%   Copyright 2011 The MathWorks, Inc.

if nargin<4 || isempty(linelabels)
    % Define blank titles if none given
    linelabels = cell(length(hLabel),1);
end

% These handles should get labeled data tips
for j = 1:length(hLabel)
    hB = hggetbehavior(hLabel(j),'datacursor');
    oldval = get(hB,'UpdateFcn');
    if isempty(oldval)
        oldval = {};
    end
    set(hB,'UpdateFcn',@(obj,evt) doLabeledDataTip(labels,linelabels{j},obj,evt,oldval{:}));
end

% But not these
if nargin>=3
    for j = 1:length(hUnlabel)
        set(hUnlabel(j),'HitTest','off');
        hB = hggetbehavior(hUnlabel(j),'datacursor');
        set(hB,'Enable',false);
    end
end

%===================== callback ====================
function datatipTxt = doLabeledDataTip(labels,linelabel,obj,evt,varargin)

ind = get(evt,'DataIndex');

% Create an observation label, or get one that already exists
if isempty(labels)
    datatipTxt = getString(message('stats:internal:addLabeledDataTip:ObservationString',num2str(ind)));
elseif isa(labels,'function_handle');
    datatipTxt = labels(ind);
else
    datatipTxt = labels{ind};
end

if isempty(varargin)
    % Put this at the top of the default data tip
    dataCursor = obj.DataCursorHandle;
    hgObject = handle(evt.Target);
    str = default_getDatatipText(dataCursor,hgObject,evt);
else
    % Put this at the top of the other data tip
    str = feval(varargin{1},obj,evt,varargin{2:end});
end

if ~isempty(str)
    if iscell(str)
        datatipTxt = [datatipTxt,sprintf('\n%s',str{:})];
    else
        datatipTxt = sprintf('%s\n%s',datatipTxt,str);
    end
end

% Add the tip title, if any
if ~isempty(linelabel)
    datatipTxt = sprintf('%s\n%s',linelabel,datatipTxt);
end
