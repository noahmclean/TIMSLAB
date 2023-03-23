classdef analyte
    %ANALYTE Superclass for analytes: reference materials and unknowns
    %   An analyte is anything you measure on a mass spectrometer
    
    properties
        name                 string   % name
        otherNames           string   % string array of alternate names
        element              string   % abbreviation for element, eg "Pb"
        nIsotopes            uint8    % number of isotopes
        species              string   % string array of species names
        relativeAbundances   double   % relative abundances, with 1 at denominator
        normalizedAbundances double   % normalized abundances (of all isotopes)

    end % properties
    
    methods

        function nIsotopes = countIsotopes(obj)
            %COUNTISOTOPES Count isotopes in analyte
            %   Only isotopes of analyte species (for now)
            nIsotopes = length(obj.relativeAbundances);
        end

        function normalizedAbundances = normalizeAbundances(obj)
            %NORMALIZEABUNDANCES Normalize abundances defined in relativeAbundances
            % sum(normalizedAbundances) = 1
            normalizedAbundances = obj.relativeAbundances/sum(obj.relativeAbundances);
        end


    end % methods

end % classdef

