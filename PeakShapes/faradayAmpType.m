classdef faradayAmpType
    %IONCOUNTERTYPE Type of amplifier for a faraday collector on a mass spectrometer
    %   resistance-based op-amp (resistance) or proprietary amp from
    %   Isotopx (ATONA)
    
    enumeration
        resistance, ATONA
    end
    
    methods
        function tf = isResistance(obj)
            %ISFARADAY Is the collector a Faraday?
            tf = obj == faradayAmpType.resistance;
        end
        
        function tf = isATONA(obj)
            %ISIONOUNTER Is the collector an ion counter?
            tf = obj == faradayAmpType.ATONA;
        end
    end % methods
    
end % classdef