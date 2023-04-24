classdef massSpecModel
    %MASSSPEC Mass Spectrometer Setup Class
    %   Properties of mass spectrometer used
    
    properties 

        collectorWidthMM = 0.85    % collector aperture width (mm)
        theoreticalBeamWidthMM  % a priori estimate of beam width (mm)
        effectiveRadiusMagnetMM % effective radius of magnet (mm)
        faradayNames            % names of Faradays as string array
        ionCounterNames         % names of ion counters as string array
        ionCounterTypes         % types of ion counters (eg, PM or EM)
        amplifierResistance     % resistance of Faraday amplifiers (ohms)
        ionCounterDeadTimes     % dead time, ns

    end

    properties (Constant)

        kB = 1.38064852e-23;   % Boltzmann constant (m^2 kg s^-2 K^-1)
        tempInK = 290;         % temperature of amplifiers (Kelvin)

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
                massSpec.theoreticalBeamWidthMM = 0.35;
                massSpec.effectiveRadiusMagnetMM = 540;
                massSpec.faradayNames = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
                massSpec.ionCounterNames = ["PM", "RS"];
                massSpec.ionCounterTypes = ["PM", "EM"];
                massSpec.ionCounterDeadTimes = [30.2, 0];
                massSpec.amplifierResistance = 1e12*ones(1,9);

                case "PhoenixKansas_1e11"
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
        

        % MORE METHODS

        function cpsPerVoltOut = cpsPerVolt(obj) % elementary charge, Ohm's Law
            cpsPerVoltOut =  6.241509074460763e+18 ./ ...
                             obj.amplifierResistance;  
        end % function cpsPerVoltOut
        
        function voltsPerCPSout = voltsPerCPS(obj) % elementary charge, Ohm's Law
            voltsPerCPSout = 1.602176634e-19 * obj.amplifierResistance; 
        end % function voltsPerCPSout

    end % methods

end % classdef