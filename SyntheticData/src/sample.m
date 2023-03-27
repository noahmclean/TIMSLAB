classdef sample < analyte
    %UNKNOWN Information about an unknown sample
    %   Unknowns are analytes that that do not have known ICs.
    
    properties
        
        % no properties yet

    end % properties
    
    methods
        function spl = sample(name, element, species, relativeAbundances)
            %UNKNOWN Construct an instance of this class
            %   Create new 'unknown' object for

            if length(species) ~= length(relativeAbundances)
                disp("BAD INPUT: # of species does not match # of abundances")
                return
            end

            spl.name = name;
            spl.element = element;
            
            % enforce letters before numbers, e.g. "Pb204" instead of "204Pb"
            spl.species = extract(species, lettersPattern) + extract(species, digitsPattern);
            
            % enforce assumption that one and only one relative abundance = 1
            % (this is the denominator/monitor/key isotope)
            if sum((relativeAbundances == 1) == 1)
                spl.relativeAbundances = relativeAbundances;
            else
                error("sampleClass:mustSpecifyDenominatorIsotope", ...
                    "One isotope in relative abundances vector must be 1.")
            end

            % calculate useful derived parameters - just use method?
            spl.nSpecies = countSpecies(spl);
            spl.normalizedAbundances = normalizeAbundances(spl);

        end % constructor function
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
        
    end % methods

end % classdef

