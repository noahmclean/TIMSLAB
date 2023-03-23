classdef unknown < analyte
    %UNKNOWN Information about an unknown sample
    %   Unknowns are analytes that that do not have known ICs.
    
    properties
        
        % no properties yet

    end % properties
    
    methods
        function unk = unknown(name, element, species, relativeAbundances)
            %UNKNOWN Construct an instance of this class
            %   Create new 'unknown' object for

            if length(species) ~= length(relativeAbundances)
                disp("BAD INPUT: # of species does not match # of abundances")
                return
            end

            unk.name = name;
            unk.element = element;
            unk.species = species;
            unk.relativeAbundances = relativeAbundances;

            % calculate useful derived parameters - just use method?
            unk.nIsotopes = countIsotopes(unk);
            unk.normalizedAbundances = normalizeAbundances(unk);

        end % constructor function
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
        
    end % methods

end % classdef

