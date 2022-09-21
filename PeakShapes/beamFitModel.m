classdef beamFitModel
    %FITSMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        G (:,:) double            % integration matrix, does convolution
        beamshape (:,1) double    % beam shape interpolated at masses beamMassInterp
        baseline (:,1) double     % baseline - when no beam in detector
        measPeakIntensityBLcorr   % measured intensity corrected for baseline
        magnetMassesWithBeam      % magnet masses for measPeakIntensityBLcorr
        measPeakIntensityWithBeam % measPeakIntensityBLcorr clipped to have beam
        leftBoundary (1,1) double % low-mass boundary of beam at thesholdIntensity
        rightBoundary (1,1) double % high-mass boundary of beam at thresholdIntensity
        GB (:,:) double           % design matrix for fit, G*B
        Wdata (:,:)               % inverse covariance matrix for data

        beamWLS   (:,1) double    % weighted least squares fit for beam parameters
        beamWNNLS (:,1) double    % nonnegative WLS fit for beam parameters
    end
    
    methods
        function beamFit = beamFitModel(data, splineBasis, peakMeas)
            %FITSMODEL Construct an instance of this class
            %   Detailed explanation goes here

            deltabeamMassInterp = splineBasis.beamMassInterp(2) - ...
                                  splineBasis.beamMassInterp(1);
            nMagnetMasses = length(data.magnetMasses);
            G = zeros(nMagnetMasses, length(splineBasis.beamMassInterp));
            for iMass = 1:nMagnetMasses % a row for each manget mass

                % massesInCollector are *model* masses
                massesInCollector = peakMeas.collectorLimits(iMass,1) <= ...
                                    splineBasis.beamMassInterp & ...
                                    splineBasis.beamMassInterp <= ...
                                    peakMeas.collectorLimits(iMass,2);

                firstMassIndexInside = find(massesInCollector,1,'first');
                lastMassIndexInside  = find(massesInCollector,1,'last');
                G(iMass, firstMassIndexInside + 1: lastMassIndexInside -1) = deltabeamMassInterp;
                G(iMass, [firstMassIndexInside, lastMassIndexInside]) = deltabeamMassInterp/2;

            end

            % calculate baseline
            hasModelBeam = any(G,2); % magnet masses with beam model mass in collector
            beamFit.baseline = data.measPeakIntensity(~hasModelBeam); % baseline when no beam in detector
            beamFit.measPeakIntensityBLcorr = data.measPeakIntensity - ...
                                       mean(data.baseline(2:end)); % 1st measured intensity is dodgy
            
            % trim data
            beamFit.G = G(hasModelBeam,:);
            beamFit.magnetMassesWithBeam = data.magnetMasses(hasModelBeam);
            beamFit.measPeakIntensityWithBeam = data.measPeakIntensityBLcorr(hasModelBeam);

            % prepare for fits - G*B (design matrix for fit) is useful to calculate upfront
            beamFit.GB = beamFit.G * splineBasis.B;

            % construct Wdata, the inverse covariance matrix for measured data
            if data.detectorName == "PhotoMultiplier"
                countData = beamFit.measPeakIntensityWithBeam * data.integPeriodMS/1000; % counts
                countData = max(countData,1)*2;  % eliminate zeros, inflate by factor 2
                beamFit.Wdata = diag(1./countData); % 
            else % if a Faraday
                beamFit.Wdata = []; % Johnson and shot noise
            end % if Photomultiplier (Daly) or Faraday

        end % constructor function
        
        function beamFit = fitLeastSquares(beamFit)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            GB = beamFit.GB; %#ok<*PROP> 
            Wdata = beamFit.Wdata;
            meas = beamFit.measPeakIntensityWithBeam;
            beamFit.beamWLS = (GB'*Wdata*GB)\(GB'*Wdata*meas);
            beamFit.beamWNNLS = lsqnonneg(chol(Wdata)*GB,chol(Wdata)*meas);

        end % function fitLeastSquares

        function beamFit = fitSmoothSpline

        end % function fitSmoothSpline

    end % methods


end % classdef

