classdef fitModel
    %FITMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        modelVector
        residualVector
        chi2
        df
        cv
        aic
    end
    
    methods
        
        function fit = fitModel(beamFitModel,inputArg2)
            %FITMODEL Construct an instance of this class
            %   Detailed explanation goes here
            fit.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end

    end % methods

end % classdef

