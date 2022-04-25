classdef massSpecModel
    %MASSSPEC Mass Spectrometer Setup Class
    %   Properties of mass spectrometer used
    
    properties

        collectorWidthMM        % collector aperture width (mm)
        theoreticalBeamWidthMM  % a priori estimate of beam width (mm)
        effectiveRadiusMagnetMM % effective radius of magnet (mm)
        faradayNames            % names of Faradays as string array
        ionCounterNames         % names of ion counters as string array
        amplifierResistance     % resistance of Faraday amplifiers (ohms)

    end
    
    methods
        function massSpec = massSpecModel(massSpecName)
            %MASSSPEC Construct an instance of this class
            %   Properties depend on input string massSpecName
            
            arguments
                massSpecName (1,1) string = ""
            end

            switch massSpecName % define properties per mass spec

                case "PhoenixKansas_1e12"
                massSpec.collectorWidthMM = 0.95135;
                massSpec.theoreticalBeamWidthMM = 0.35;
                massSpec.effectiveRadiusMagnetMM = 540;
                massSpec.faradayNames = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
                massSpec.ionCounterNames = ["PM", "SEM"];
                massSpec.amplifierResistance = 1e12;

                case "PhoenixKansas_1e11"
                massSpec.collectorWidthMM = 0.95135;
                massSpec.theoreticalBeamWidthMM = 0.35;
                massSpec.effectiveRadiusMagnetMM = 540;
                massSpec.faradayNames = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
                massSpec.ionCounterNames = ["PM", "SEM"];
                massSpec.amplifierResistance = 1e11;

                case ""
                disp(" ")
                disp("Mass spectrometer not recognized")
                return

            end % switch massSpecName
        
        end % constructor function
        
        function collectorWidthAMU = calcCollectorWidthAMU(massSpec, massAtCenterAMU)
            %CALCCOLLECTORWIDTHAMU Collector width in AMU
            %   Calculate collector aperture width in AMU at measured mass
            %   massAtCenterAMU is average/peak mass of the scan, in amu
            
            collectorWidthAMU = massAtCenterAMU / ...
                massSpec.effectiveRadiusMagnetMM * massSpec.collectorWidthMM;

        end % calcCollectorWidthAMU

        function beamWidthAMU = calcBeamWidthAMU(massSpec, massAtCenterAMU)
            %CALCBEAMWIDTHAMU Beam width in AMU
            %   Calculate beam width in AMU at measured mass
            %   massAtCenterAMU is average/peak mass of the scan, in amu

            beamWidthAMU = massAtCenterAMU / ...
                massSpec.effectiveRadiusMagnetMM * massSpec.theoreticalBeamWidthMM;
        end % calcBeamWidthAMU

    end % methods

end % classdef