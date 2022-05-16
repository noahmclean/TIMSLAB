classdef peakMeasProperties
    %PEAKMEASPROPERTIES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        collectorWidthAMU       % collector width in AMU
        theoreticalBeamWidthAMU % theoretical beam width in AMU
        collectorLimits         % mass range in collector at each magnet mass 
        deltaMagnetMass         % change in magnet mass between measurements
        beamWindow              % [min, max] mass range to model beam over
    end
    
    methods
        function peakMeas = peakMeasProperties(data,massSpec)

            %PEAKMEASPROPERTIES Construct an instance of this class
            %   Calculations about the peak measurement on this mass spec
            
            if nargin > 0 % if called with arguments

            peakMeas.collectorWidthAMU       = calcWidthInAMU(data, massSpec, massSpec.collectorWidthMM);
            peakMeas.theoreticalBeamWidthAMU = calcWidthInAMU(data, massSpec, massSpec.theoreticalBeamWidthMM);
            
            % collectorLimits is a matrix with two columns and the same
            % number of rows as magnet masses.  Each row contains the mass
            % range of the beam that is entering the collector (defined by
            % collectorWidthAMU)
            peakMeas.collectorLimits = data.magnetMasses + ...
                      [-data.collectorWidthAMU, data.collectorWidthAMU]/2;
            
            
            peakMeas.deltaMagnetMass = data.magnetMasses(2)-data.magnetMasses(1);
            
            peakMeas.beamWindow = data.theoreticalBeamWidthAMU*2;

            end

        end % constructor function
        
    end % public methods block

end % classdef

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function widthAMU = calcWidthInAMU(data, massSpec, widthMM)

widthAMU = data.peakCenterMass / ...
          massSpec.effectiveRadiusMagnetMM * ...
          widthMM;

end % function