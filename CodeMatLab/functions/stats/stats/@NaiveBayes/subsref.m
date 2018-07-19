function [varargout] = subsref(obj,s)
%SUBSREF Subscripted reference for a NaiveBayes object.
%   B = SUBSREF(OBJ,S) is called for the syntax OBJ(S) when OBJ is a
%   NaiveBayes object. S is a structure array with the fields:
%       type -- string containing '()', '{}', or '.' specifying the
%               subscript type.
%       subs -- Cell array or string containing the actual subscripts.
%
%   See also NAIVEBAYES.

%   Copyright 2008 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2010/10/08 17:27:49 $

switch s(1).type
    case '()'
        error(message('stats:NaiveBayes:subsref:ArraySubscript', class( obj )));
        
    case '.'
        methodsProp=[methods(obj);properties(obj)];
        if ~any(strcmp(s(1).subs, methodsProp))
            error(message('stats:NaiveBayes:subsref:AccessPrivate', s( 1 ).subs, class( obj )));
        elseif   strcmp(s(1).subs,'fit')
            error(message('stats:NaiveBayes:subsref:InstCallStatic', s( 1 ).subs, class( obj )));
        end
        if isequal(s(1).subs,'display') %% nargout==0
            display(obj,inputname(1));
        else
            [varargout{1:nargout}] = builtin('subsref',obj,s);
        end
    otherwise
        [varargout{1:nargout}] = builtin('subsref',obj,s);
end
