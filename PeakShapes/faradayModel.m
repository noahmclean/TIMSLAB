classdef faradayModel < collectorModel
    %IONCOUNTERMODEL an ion counter
    %   Detailed explanation goes here
    
    properties
        ampType                faradayAmpType   % resistance or ATONA
        amplifierResistance    (1,1) double     % unit: ohms
    end
    
    methods
        function obj = faradayModel(code, name, ampType, amplifierResistance)
            %IONCOUNTERMODEL Construct an instance of this class
            %   An ionCounterModel goes into a collectorArray in a
            %   massSpecModel
            obj@collectorModel(code, name)
            obj.ampType = ampType;
            obj.amplifierResistance = amplifierResistance;
        end
        
        function tf = isATONA(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            tf = obj.ampType == faradayAmpType.ATONA;
        end

        function tf = isResistance(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            tf = obj.ampType == faradayAmpType.resistance;
        end

    end
end

