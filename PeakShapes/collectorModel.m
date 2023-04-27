classdef collectorModel < matlab.mixin.Heterogeneous
    %COLLECTORModel Properties of a mass spectrometer collector
    %   Superclass for ion counters and faradays
    
    properties
        code string %e.g., "H2"    or "PM"
        name string %e.g., "High2" or "Photomultiplier"
    end
    
    methods
        function obj = collectorModel(code, name)
            %COLLECTORMODEL constructor method given code and name
            obj.code = code;
            obj.name = name;
        end
        
    end % methods

end % classdef