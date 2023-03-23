classdef referenceMaterial < analyte
    %REFERENCEMATERIAL Relevant Reference Material Data
    %   Isotopic Composition Reference Materials
    
    properties

        reference            string   % reference

    end
    
    methods
        function refmat = referenceMaterial(name)
            %REFERENCEMATERIAL Construct an instance of this class
            %   refmat 
            refmat.name = name;

        switch name

            case "NBS981"
                refmat.otherNames         = ["NBS 981", "981"];
                refmat.element            = "Pb";
                refmat.species            = ["Pb204", "Pb206", "Pb207", "Pb208"];
                refmat.relativeAbundances = [0.0590074 1       0.914683 2.1681];
                refmat.reference          = "Condon et al. (2015)";

            case "NBS982"
                refmat.otherNames         = ["NBS 982", "982"];
                refmat.element            = "Pb";
                refmat.species            = ["Pb204",  "Pb206", "Pb207", "Pb208"];
                refmat.relativeAbundances = [0.0272058 1        0.466967 1.000249];
                refmat.reference          = "Condon et al. (2015)";

            case "SRM987"
                refmat.otherNames         = ["NBS987", "SRM 987", "NBS 987"];
                refmat.element            = "Sr";
                refmat.species            = ["Sr84",  "Sr86", "Sr87",  "Sr88"];
                refmat.relativeAbundances = [0.056493 1       0.710249 1/0.119400];
                refmat.reference          = "Brennecka et al. (2013)";

            case "SmKU"
                refmat.otherNames         = ["Ames Sm", "KU Sm", "KUSm", "Sm KU", "SmKU1"];
                refmat.element            = "Sm";
                refmat.species            = ["Sm144", "Sm147", "Sm148", "Sm149", "Sm150", "Sm152", "Sm154"];
                refmat.relativeAbundances = [0.416583 2.031957 1.523370 1.872696 1        3.623281 3.082655];
                refmat.reference          = "Brennecka et al. (2013)";

            case "JNdi"
                refmat.otherNames         = "JNdi Nd";
                refmat.element            = "Nd";
                refmat.species            = ["Nd142", "Nd143", "Nd144", "Nd145", "Nd146", "Nd148", "Nd150"];
                refmat.relativeAbundances = [1.141834 0.512098 1        0.348402 0.7219   0.241579 0.236455];
                refmat.reference          = "Brennecka et al. (2013)";

        end % switch case name

        % calculate useful derived parameters - just use method?
        refmat.nSpecies = countSpecies(refmat);
        refmat.normalizedAbundances = normalizeAbundances(refmat);

        end % constructor function
        
    end % methods
end % classdef

