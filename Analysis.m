%to be able to create an array of Analysis, 
%and then assign K0SAnalysis and LambdaAnalysis objects to array elements
classdef Analysis < handle & matlab.mixin.Heterogeneous
    properties

    end
    methods
        function obj = Analysis()
        end
        function start(obj)
        end
        function stop(obj)
        end
        function event(obj)
        end
    end
end