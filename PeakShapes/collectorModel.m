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
    end

    methods(Sealed) % can't be overridden by a subclass

        function tf = isFaraday(coll)
            nCollectors = length(coll);
            tf = false(nCollectors,1);
            for iColl = 1:nCollectors
                tf(iColl) = isa(coll(iColl), 'faradayModel');
            end %for
        end % isFaraday

        function tf = isIonCounter(coll)
            nCollectors = length(coll);
            tf = false(nCollectors,1);
            for iColl = 1:nCollectors
                tf(iColl) = isa(coll(iColl), 'ionCounterModel');
            end %for
        end % isIonCounter
        
    end % methods(Sealed)

end % classdef