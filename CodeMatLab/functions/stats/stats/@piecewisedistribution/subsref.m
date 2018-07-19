function varargout = subsref(t,s)
%SUBSREF Subscripted reference for a piecewisedistribution object.
%   B = SUBSREF(T,S) is called for the syntax T(X) when T is a
%   piecewisedistribution object.  S is a structure array with the fields:
%       type -- string containing '()', '{}', or '.' specifying the
%               subscript type.
%       subs -- Cell array or string containing the actual subscripts.

%   Copyright 2006-2007 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2011/05/09 01:28:18 $

switch s(1).type
case '()'
    if ~isscalar(s)
        error(message('stats:piecewisedistribution:subsref:InvalidSubscriptExpr'));
    
    elseif numel(s(1).subs) ~= 1
        error(message('stats:piecewisedistribution:subsref:ParenArgNotSingleArray'));
    end
    
    % For now we allow this form to invoke the cdf method, but warn that
    % this may change in the future
    warning(message('stats:piecewisedistribution:subsref:ObsoleteSyntax'));

    % '()' is a reference to the cdf method
    [varargout{1:nargout}] = cdf(t,s(1).subs{1});


case '{}'
    error(message('stats:piecewisedistribution:subsref:CellSubscript'));
        
case '.'
    % Support dot subscripting to call methods
    if isscalar(s)
        args = {};
    elseif numel(s)>2 || ~isequal(s(2).type,'()')
         error(message('stats:piecewisedistribution:subsref:DotSubscriptWithNonScalarsNotAllowed'));
    else   % numel(s)==2
        args = s(2).subs;
    end
    
    % For improved performances these methods are hard-coded instead of
    % being found via the "methods" function and called via "feval"
    switch(s(1).subs)
      case 'cdf'
        [varargout{1:nargout}] = cdf(t,args{:});
      case 'icdf'
        [varargout{1:nargout}] = icdf(t,args{:});
      case 'random'
        [varargout{1:nargout}] = random(t,args{:});
      case 'pdf'
        [varargout{1:nargout}] = pdf(t,args{:});
      case 'segment'
        [varargout{1:nargout}] = segment(t,args{:});
      case 'nsegments'
        [varargout{1:nargout}] = nsegments(t,args{:});
      case 'boundary'
        [varargout{1:nargout}] = boundary(t,args{:});
      case 'disp'
        [varargout{1:nargout}] = disp(t,args{:});
        
      otherwise
         error(message('stats:piecewisedistribution:subsref:BadMethod', s( 1 ).subs));
    end
end
