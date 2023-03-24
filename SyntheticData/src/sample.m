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
            spl.species = extract(species, lettersPattern) + extract(species, digitsPattern);
            spl.relativeAbundances = relativeAbundances;

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

