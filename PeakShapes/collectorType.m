classdef collectorType
    %COLLECTORTYPE Type of mass spectrometer collector
    %   collectors are ion counters or faradays
    
    enumeration
        faraday, ionCounter
    end
    
    methods
        function tf = isFaraday(obj)
            %ISFARADAY Is the collector a Faraday?
            tf = obj == collectorType.faraday;
        end
        
        function tf = isIonCounter(obj)
            %ISIONOUNTER Is the collector an ion counter?
            tf = obj == collectorType.ionCounter;
        end
    end % methods

end % classdef