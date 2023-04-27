classdef faradayModel < collectorModel
    %IONCOUNTERMODEL an ion counter
    %   Detailed explanation goes here
    
    properties
        ampType                faradayAmpType
        amplifierResistance    (1,1) double
    end
    
    methods
        function obj = faradayModel(code, name)
            %IONCOUNTERMODEL Construct an instance of this class
            %   Detailed explanation goes here
            obj@collectorModel(code, name)
        end
        
        function tf = isATONA(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            tf = obj.type + faradayAmpType.ATONA;
        end
    end
end

