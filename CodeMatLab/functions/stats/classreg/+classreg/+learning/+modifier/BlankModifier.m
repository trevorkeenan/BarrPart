classdef BlankModifier < classreg.learning.modifier.Modifier

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2011/09/19 18:11:56 $
    
    properties(Constant=true,GetAccess=public)
        FitInfoDescription = getString(message('stats:classreg:learning:modifier:BlankModifier:FitInfoDescription'));
    end
    
    methods
        function this = BlankModifier()
            this = this@classreg.learning.modifier.Modifier(0,1);
        end
        
        function [this,mustTerminate,X,Y,W,fitData] = modify(this,X,Y,W,H,fitData)
            mustTerminate = false;
        end
        
        function combiner = makeCombiner(this)
            combiner = classreg.learning.combiner.WeightedAverage(ones(this.T,1));
        end
    end
    
end
