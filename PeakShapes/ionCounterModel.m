classdef ionCounterModel < collectorModel
    %IONCOUNTERMODEL an ion counter
    %   Detailed explanation goes here
    
    properties
        type        ionCounterType
        deadtime    (1,1) double
        darknoise   (1,1) double
    end
    
    methods
        function obj = ionCounterModel(code, name)
            %IONCOUNTERMODEL Construct an instance of this class
            %   Detailed explanation goes here
            obj@collectorModel(code, name)
        end
        
        function tf = isDaly(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            tf = obj.type + ionCounterType.Daly;
        end
    end
end

