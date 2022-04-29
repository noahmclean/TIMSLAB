classdef peakMeasProperties
    %PEAKMEASPROPERTIES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        collectorWidthAMU       % collector width in AMU
        theoreticalBeamWidthAMU % theoretical beam width in AMU
        collectorLimits         % mass range in collector at each magnet mass 
        deltaMangetMass         % change in magnet mass between measurements
        beamWindow              % [min, max] mass range to model beam over
    end
    
    methods
        function peakMeas = peakMeasProperties(data,massSpec)
            %PEAKMEASPROPERTIES Construct an instance of this class
            %   Calculations about the peak measurement on this mass spec
            
            peakMeas.collectorWidthAMU       = calcWidthInAMU(data, massSpec, massSpec.collectorWidthMM);
            peakMeas.theoreticalBeamWidthAMU = calcWidthInAMU(data, massSpec, massSpec.theoreticalBeamWidthMM);
            
            

        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function widthAMU = calcWidthInAMU(data, massSpec, widthMM)

widthAMU = data.peakCenterMass / ...
          massSpec.effectiveRadiusMagnetMM * ...
          widthMM;

end % function