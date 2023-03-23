classdef massSpecModel
    %MASSSPEC Mass Spectrometer Setup Class
    %   Properties of mass spectrometer used
    
    properties 

        collectorWidthMM        % collector aperture width (mm)
        theoreticalBeamWidthMM  % a priori estimate of beam width (mm)
        effectiveRadiusMagnetMM % effective radius of magnet (mm)
        faradayNames            % names of Faradays as string array
        ionCounterNames         % names of ion counters as string array
        ionCounterTypes         % types of ion counters (eg, PM or EM)
        amplifierResistance     % resistance of Faraday amplifiers (ohms)
        ionCounterDeadTimes     % dead time, ns

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
                %massSpec.collectorWidthMM = 0.95135;
                massSpec.collectorWidthMM = 0.855;
                massSpec.theoreticalBeamWidthMM = 0.35;
                massSpec.effectiveRadiusMagnetMM = 540;
                massSpec.faradayNames = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
                massSpec.ionCounterNames = ["PM", "RS"];
                massSpec.ionCounterTypes = ["PM", "EM"];
                massSpec.ionCounterDeadTimes = [30.2, 0];
                massSpec.amplifierResistance = 1e12*ones(1,9);

                case "PhoenixKansas_1e11"
                massSpec.collectorWidthMM = 0.95135;
                massSpec.theoreticalBeamWidthMM = 0.35;
                massSpec.effectiveRadiusMagnetMM = 540;
                massSpec.faradayNames = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
                massSpec.ionCounterNames = ["PM", "RS"];
                massSpec.ionCounterTypes = ["PM", "EM"];
                massSpec.ionCounterDeadTimes = [30.2, 0];
                massSpec.amplifierResistance = 1e11*ones(1,9);

                case ""
                disp(" ")
                disp("Mass spectrometer not recognized")
                return

            end % switch massSpecName
        
        end % constructor function
        


    end % methods

end % classdef