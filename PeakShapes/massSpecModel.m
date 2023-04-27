classdef massSpecModel
    %MASSSPEC Mass Spectrometer Setup Class
    %   Properties of mass spectrometer used
    
    properties 

        nominalCollectorWidthMM = 1 % collector aperture width (mm)
        theoreticalBeamWidthMM  % a priori estimate of beam width (mm)
        effectiveRadiusMagnetMM % effective radius of magnet (mm)
        faradayNames            % names of Faradays as string array, eg "High2"
        faradayCodes            % codes for Faradays, eg "H2"
        faradayTypes            % "resistance" or "ATONA"
        faradayApertureWidthMM  % aperture widths for Faraday collectors (~1 mm)
        amplifierResistance     % resistance of Faraday amplifiers (ohms)
        ionCounterNames         % names of ion counters as string array
        ionCounterTypes         % types of ion counters (eg, PM or EM)
        ionCounterDeadTimes     % dead time, ns
        ionCounterDarkNoise     % dark noise, cps
        ionCounterAperturewidthMM % aperture widths for ion counters (~1 mm)

    end

    properties (Constant)

        kB = 1.38064852e-23;   % Boltzmann constant (m^2 kg s^-2 K^-1)
        tempInK = 290;         % temperature of amplifiers (Kelvin)
        coulomb = 6241509074460762607.776; % 1 coulomb in elementary charges
        % Ickert: "1/t * 1484 = amplifier noise, t in seconds, noise in cps"
        noiseConstantATONAsCPS = 1484; % cps noise

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
                massSpec.faradayNames = ["Low5", "Low4", "Low3", "Low2", "Axial", "High1", "High2", "High3", "High4"];
                massSpec.faradayCodes = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
                massSpec.faradayTypes = repmat("resistance", 1, 9);
                massSpec.amplifierResistance = 1e12*ones(1,9);
                massSpec.ionCounterNames = ["PM", "RS"];
                massSpec.ionCounterTypes = ["PM", "EM"];
                massSpec.ionCounterDeadTimes = [30.2, 0];
                massSpec.ionCounterDarkNoise = 0.1;
                

                case "PhoenixKansas_1e11"
                massSpec.theoreticalBeamWidthMM = 0.35;
                massSpec.effectiveRadiusMagnetMM = 540;
                massSpec.faradayNames = ["Low5", "Low4", "Low3", "Low2", "Axial", "High1", "High2", "High3", "High4"];
                massSpec.faradayCodes = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
                massSpec.faradayTypes = repmat("resistance", 1, 9);
                massSpec.amplifierResistance = 1e11*ones(1,9);
                massSpec.ionCounterNames = ["PM", "RS"];
                massSpec.ionCounterTypes = ["PM", "EM"];
                massSpec.ionCounterDeadTimes = [30.2, 0];
                massSpec.ionCounterDarkNoise = 0.1;
                

                case "PhoenixPurdue_ATONA"
                massSpec.theoreticalBeamWidthMM = 0.35;
                massSpec.effectiveRadiusMagnetMM = 540;
                massSpec.faradayNames = ["Low5", "Low4", "Low3", "Low2", "Axial", "High1", "High2", "High3", "High4"];
                massSpec.faradayCodes = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
                massSpec.faradayTypes = repmat("ATONA", 1, 9);
                massSpec.amplifierResistance = 1e11*ones(1,9); % results given in volts of 10^11-equivalent
                massSpec.ionCounterNames = "PM";
                massSpec.ionCounterTypes = "PM";
                massSpec.ionCounterDeadTimes = [31.3, 0];
                massSpec.ionCounterDarkNoise = 0.1;

                case ""
                disp(" ")
                disp("Mass spectrometer not recognized")
                return

            end % switch massSpecName
        
        end % constructor function
        

        % MORE METHODS

        function cpsPerVoltOut = cpsPerVolt(obj) % elementary charge, Ohm's Law
            cpsPerVoltOut =  obj.coulomb ./ ...
                             obj.amplifierResistance;  
        end % function cpsPerVoltOut
        
        function voltsPerCPSout = voltsPerCPS(obj) % elementary charge, Ohm's Law
            voltsPerCPSout =  obj.amplifierResistance ./ ...
                              obj.coulomb; 
        end % function voltsPerCPSout

    end % methods

end % classdef