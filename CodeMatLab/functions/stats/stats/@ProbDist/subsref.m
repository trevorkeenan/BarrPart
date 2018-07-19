function [varargout] = subsref(a,s)
%SUBSREF Subscripted reference for a probability distribution.
%   B = SUBSREF(A,S) is called for the syntax A.PROP where A is a
%   parametric probability distribution and PROP is a property name (such
%   as a parameter name).   S is a structure array with the fields:
%       type -- string containing '.' specifying the subscript type.
%       subs -- Cell array or string containing the parameter name.
%
%   B = A.PROP or A.(PROPNAME) returns the value of the specified
%   property for the distribution A.  PROP is a property name literal,
%   or PROPNAME is a character variable containing a property name.
%
%   B = A(I,J) and B = A{I,J} are not allowed.
%
%   See also ProbDist, ProbDistUnivParam.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2011/05/09 01:27:45 $

switch s(1).type
case '()'
    if isscalar(s) && isequal(s(1).subs,{1 1})
        varargout{1} = a;
    else
        error(message('stats:ProbDist:subsref:NotAllowed'));
    end
case '{}'
    error(message('stats:ProbDist:subsref:CellNotAllowed'));
        
case '.'
    % A reference to a property name or parameter name.  Could be any sort
    % of subscript following that.
    pname = s(1).subs;

    if ischar(pname) && size(pname,1) == 1
        if any(strcmp(pname,properties(a)))
            % Property name
            b = a.(pname);
        elseif any(strcmp(pname,methods(a)))
            % Method call via dot subscripting
            if isscalar(s)
                args = {};
            elseif numel(s)>2 || ~isequal(s(2).type,'()')
                error(message('stats:ProbDist:subsref:BadSubscript'));
            else   % numel(s)==2
                args = s(2).subs;
            end
            [varargout{1:nargout}] = feval(s(1).subs,a,args{:});
            return
        else
            error(message('stats:ProbDist:subsref:BadProperty', pname));
        end
    else
        error(message('stats:ProbDist:subsref:IllegalSubscript'));
    end

    if ~isscalar(s)
        % Handle cascaded subscripting
        [varargout{1:nargout}] = subsref(b,s(2:end));
    else
        % No cascaded subscripting
        varargout{1} = b;
    end
end
