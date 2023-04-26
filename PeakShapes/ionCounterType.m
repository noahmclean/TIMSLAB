classdef ionCounterType
    %IONCOUNTERTYPE Type of ion counter on a mass spectrometer
    %   ion counters are Dalys or SEMs (secondary electron multipliers)
    
    enumeration
        Daly, SEM
    end
    
    methods
        function tf = isDaly(obj)
            %ISFARADAY Is the collector a Faraday?
            tf = obj == ionCounterType.Daly;
        end
        
        function tf = isSEM(obj)
            %ISIONOUNTER Is the collector an ion counter?
            tf = obj == ionCounterType.SEM;
        end
    end % methods
    
end % classdef