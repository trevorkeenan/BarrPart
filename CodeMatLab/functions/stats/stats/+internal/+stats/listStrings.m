function s = liststrings(c)
% LISTSTRINGS Create a list from a cell array of strings.
%    T = LISTSTRINGS(C) takes the strings in the cell array C and creates a
%    comma-separated list from them.

%   $Revision: 1.1.6.2 $  $Date: 2011/09/19 18:13:08 $
%   Copyright 2011 The MathWorks, Inc.

s = sprintf(', %s',c{:});
s(1:2) = [];  % remove leading comma
