classdef ionCounterModel < collectorModel
    %IONCOUNTERMODEL an ion counter
    %   Detailed explanation goes here
    
    properties
        type        ionCounterType   % Daly or SEM
        deadtime    (1,1) double     % dead time, nanoseconds
        darknoise   (1,1) double     % dark noise, cps
    end
    
    methods
        function obj = ionCounterModel(code, name, type, deadtime, darknoise)
            %IONCOUNTERMODEL Construct an instance of this class
            %   An ionCounterModel goes into a collectorArray table in a
            %   massSpecModel 
            obj@collectorModel(code, name)
            obj.type = type;
            obj.deadtime = deadtime;
            obj.darknoise = darknoise;
        end
        
        function tf = isDaly(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            tf = obj.type == ionCounterType.Daly;
        end

        function tf = isSEM(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            tf = obj.type == ionCounterType.SEM;
        end

    end
end

