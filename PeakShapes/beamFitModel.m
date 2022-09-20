classdef beamFitModel
    %FITSMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        G (:,:) double           % integration matrix, does convolution
        beamshape {mustBeVector} % beam shape interpolated at masses beamMassInterp
        leftBoundary (1,1)       % low-mass boundary of beam at thesholdIntensity
        rightBoundary (1,1)      % high-mass boundary of beam at thresholdIntensity


        

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

