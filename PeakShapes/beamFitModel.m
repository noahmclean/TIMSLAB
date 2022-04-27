classdef beamFitModel
    %FITSMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        G (:,:) double

    end
    
    methods
        function fits = beamFitModel(data,splineBasis)
            %FITSMODEL Construct an instance of this class
            %   Detailed explanation goes here
            
        

        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

