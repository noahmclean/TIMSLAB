classdef massSpecModel
    %MASSSPEC Mass Spectrometer Setup Class
    %   Properties of mass spectrometer used
    
    properties 

        collectorArray          table   % table with collector info        
        theoreticalBeamWidthMM  double  % a priori estimate of beam width (mm)
        effectiveRadiusMagnetMM double  % effective radius of magnet (mm)

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
                PM  = ionCounterModel("PM", "Photomultiplier", "Daly", 30.1, 0.1);
                SEM = ionCounterModel("SEM", "Rear SEM", "SEM", NaN, NaN);
                L5  = faradayModel("L5", "Low5" , "resistance", 10^12);
                L4  = faradayModel("L4", "Low4" , "resistance", 10^12);
                L3  = faradayModel("L3", "Low3" , "resistance", 10^12);
                L2  = faradayModel("L2", "Low2" , "resistance", 10^12);
                Ax  = faradayModel("Ax", "Axial", "resistance", 10^12);
                H1  = faradayModel("H1", "High1", "resistance", 10^12);
                H2  = faradayModel("H2", "High2", "resistance", 10^12);
                H3  = faradayModel("H3", "High3", "resistance", 10^12);
                H4  = faradayModel("H4", "High4", "resistance", 10^12);
                collectors = [PM; SEM; L5; L4; L3; L2; Ax; H1; H2; H3; H4];
                apertureWidthMMnominal = ones(11,1);
                apertureWidthMMmeas    = nan(11,1);
                massSpec.collectorArray = table(collectors, ...
                               apertureWidthMMnominal, apertureWidthMMmeas, ...
                               'RowNames', ["PM"; "SEM"; "L5"; "L4"; "L3"; ...
                                            "L2"; "Ax"; "H1"; "H2"; "H3"; "H4"]);
                massSpec.theoreticalBeamWidthMM = 0.35;
                massSpec.effectiveRadiusMagnetMM = 540;

                case "PhoenixKansas_1e11"
                PM  = ionCounterModel("PM", "Photomultiplier", "Daly", 30.1, 0.1);
                SEM = ionCounterModel("SEM", "Rear SEM", "SEM", NaN, NaN);
                L5  = faradayModel("L5", "Low5" , "resistance", 10^11);
                L4  = faradayModel("L4", "Low4" , "resistance", 10^11);
                L3  = faradayModel("L3", "Low3" , "resistance", 10^11);
                L2  = faradayModel("L2", "Low2" , "resistance", 10^11);
                Ax  = faradayModel("Ax", "Axial", "resistance", 10^11);
                H1  = faradayModel("H1", "High1", "resistance", 10^11);
                H2  = faradayModel("H2", "High2", "resistance", 10^11);
                H3  = faradayModel("H3", "High3", "resistance", 10^11);
                H4  = faradayModel("H4", "High4", "resistance", 10^11);
                collectors = [PM; SEM; L5; L4; L3; L2; Ax; H1; H2; H3; H4];
                apertureWidthMMnominal = ones(11,1);
                apertureWidthMMmeas    = nan(11,1);
                massSpec.collectorArray = table(collectors, ...
                               apertureWidthMMnominal, apertureWidthMMmeas, ...
                               'RowNames', ["PM"; "SEM"; "L5"; "L4"; "L3"; ...
                                            "L2"; "Ax"; "H1"; "H2"; "H3"; "H4"]);
                massSpec.theoreticalBeamWidthMM = 0.35;
                massSpec.effectiveRadiusMagnetMM = 540;
                

                case "PhoenixPurdue_ATONA"
                PM  = ionCounterModel("PM", "Photomultiplier", "Daly", 30.1, 0.1);
                L5  = faradayModel("L5", "Low5" , "ATONA", 10^11);
                L4  = faradayModel("L4", "Low4" , "ATONA", 10^11);
                L3  = faradayModel("L3", "Low3" , "ATONA", 10^11);
                L2  = faradayModel("L2", "Low2" , "ATONA", 10^11);
                Ax  = faradayModel("Ax", "Axial", "ATONA", 10^11);
                H1  = faradayModel("H1", "High1", "ATONA", 10^11);
                H2  = faradayModel("H2", "High2", "ATONA", 10^11);
                H3  = faradayModel("H3", "High3", "ATONA", 10^11);
                H4  = faradayModel("H4", "High4", "ATONA", 10^11);
                collectors = [PM; L5; L4; L3; L2; Ax; H1; H2; H3; H4];
                nCollectors = length(collectors);
                apertureWidthMMnominal = ones(nCollectors,1);
                apertureWidthMMmeas    = nan(nCollectors,1);
                massSpec.collectorArray = table(collectors, ...
                               apertureWidthMMnominal, apertureWidthMMmeas, ...
                               'RowNames', ["PM"; "L5"; "L4"; "L3"; ...
                                            "L2"; "Ax"; "H1"; "H2"; "H3"; "H4"]);
                massSpec.theoreticalBeamWidthMM = 0.35;
                massSpec.effectiveRadiusMagnetMM = 540;

                otherwise
                errID = "MassSpecModel:MassSpecNameNotRecognized";
                msg = "Mass Spectrometer model name not recognized.";
                throw(MException(errID, msg))

            end % switch massSpecName
        
        end % constructor function
        

        % MORE METHODS

        % cpsPerVolt and voltsPerCPS use elementary charge, Ohm's Law

        function cpsPerVolt = cpsPerVolt(obj, collectorCode) % elementary charge, Ohm's Law
            
            idx = getCollectorIndex(obj, collectorCode);
            % get the amplifier resistance (real or voltage-equivalent)
            amplifierResistance = obj.collectorArray.collectors(idx).amplifierResistance;
            cpsPerVolt =  obj.coulomb ./ ...
                             amplifierResistance;  
        end % function cpsPerVoltOut
        
        function voltsPerCPS = voltsPerCPS(obj, collectorCode) 
            
            idx = getCollectorIndex(obj, collectorCode);
            % get the amplifier resistance (real or voltage-equivalent)
            amplifierResistance = obj.collectorArray.collectors(idx).amplifierResistance;

            voltsPerCPS =  amplifierResistance ./ ...
                              obj.coulomb; 
        end % function voltsPerCPSout

        function awMM = getApertureWidthMM(obj, codename)
            
            % convert collector code or name to table index
            idx = getCollectorIndex(obj, codename);

            % if there is no measured aperture width, use the nominal width
            if isnan(obj.collectorArray{idx, "apertureWidthMMmeas"})
                awMM = obj.collectorArray{idx, "apertureWidthMMnominal"};
            % else, use the measured aperture width
            else 
                awMM = obj.collectorArray{idx, "apertureWidthMMmeas"};
            end % if
            
        end % function get.apertureWidthMM

    end % public methods

    methods (Access = private)
        
        % return collector index given collector code or name
        function idx = getCollectorIndex(obj, codename)

            idx = [find([obj.collectorArray.collectors.code] == codename);
                   find([obj.collectorArray.collectors.name] == codename)];
            
        end

    end % private methods

end % classdef