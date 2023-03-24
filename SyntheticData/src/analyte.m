classdef analyte
    %ANALYTE Superclass for analytes: reference materials and unknowns
    %   An analyte is anything you measure on a mass spectrometer
    
    properties
        name                 string   % name
        otherNames           string   % string array of alternate names
        element              string   % abbreviation for element, eg "Pb"
        species              string   % string array of species names
        nSpecies             uint8    % number of species present
        relativeAbundances   double   % relative abundances, with 1 at denominator
        normalizedAbundances double   % normalized abundances (of all isotopes)

    end % properties
    
    methods

        function nSpecies = countSpecies(obj)
            %COUNTISOTOPES Count isotopes in analyte
            %   Only isotopes of analyte species (for now)
            nSpecies = length(obj.species);
        end

        function normalizedAbundances = normalizeAbundances(obj)
            %NORMALIZEABUNDANCES Normalize abundances defined in relativeAbundances
            % sum(normalizedAbundances) = 1
            normalizedAbundances = obj.relativeAbundances/sum(obj.relativeAbundances);
        end

        function logRatios = logRatioAbundances(obj)
            logRatios = log(obj.relativeAbundances);
        end

        function denIsoIndex = denominatorIsotopeIndex(obj)
            %DENOMINATORISOTOPEINDEX determine index of denominator isotope
            % e.g. 206 in Pb IC [204/206 207/206 208/206] 
            % assumed to be the abundance = 1 in relativeAbundances
            denIsoIndex = find(obj.relativeAbundances == 1);
        end

        function massVec = massVector(obj)
            %MASSVECTOR vector of isotopic masses of all isotopes in "species"
            % uses mass class in collectorRelativeEfficiencies
            nSpec = obj.countSpecies;
            massVec = zeros(1, nSpec);
            for iSpecies = 1:nSpec
                massVec(iSpecies) = mass.(obj.species(iSpecies));
            end
        end

        function normMassVector = normMasses(obj)
            %NORMMASSES vector of isotopic masses normalized by denominator species
            % uses mass class in collectorRelativeEfficiencies
            denIsoIndex = denominatorIsotopeIndex(obj);
            massVec = massVector(obj);
            denMass = massVec(denIsoIndex);
            normMassVector = massVec / denMass;
        end

        function logNormMassVector = logNormMasses(obj)
            %LOGNORMMASSES vector of log-ratios of masses
            % with denominatorIsotope in denominator (=1 in relative abundance vector)
            logNormMassVector = log(obj.normMasses);
        end


    end % methods

end % classdef

