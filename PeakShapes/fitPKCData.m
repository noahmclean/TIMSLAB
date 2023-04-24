function resnorm = fitPKCData(massSpec, data, collectorWidthMM)
% FITPKCDATA fit a beam shape to a peakcenter file 
% output the misfit for a given collector width
% use this function to solve for optimum collectorWidthMM

% to optimize: collector width (mm)
massSpec.collectorWidthMM = collectorWidthMM;

data.collectorWidthAMU = calcCollectorWidthAMU(data, massSpec);
%data.collectorWidthAMU = 0.139;
data.theoreticalBeamWidthAMU = calcBeamWidthAMU(data, massSpec);
peakMeas = peakMeasProperties(data, massSpec);

% spline basis B
bdeg = 3; % order of spline (= order of polynomial pieces)
pord = 2; % order of differences for pSplines
beamKnots = ceil(peakMeas.beamWindow/(peakMeas.deltaMagnetMass)) - 2*bdeg; % 3 xtra knots for cubic spline
nInterp = 1000; % number of interpolated segments

xl = data.peakCenterMass - peakMeas.beamWindow/2;
xr = data.peakCenterMass + peakMeas.beamWindow/2;
beamMassInterp = linspace(xl, xr, nInterp);

splineBasis = splineBasisModel(beamMassInterp,beamKnots,bdeg);
deltabeamMassInterp = beamMassInterp(2)-beamMassInterp(1);

% calculate integration matrix G, depends on B, data
nMagnetMasses = length(data.magnetMasses);
G = zeros(nMagnetMasses, nInterp);
for iMass = 1:nMagnetMasses % a row for each manget mass

    % massesInCollector are *model* masses
    massesInCollector = peakMeas.collectorLimits(iMass,1) <= beamMassInterp & ...
                        beamMassInterp <= peakMeas.collectorLimits(iMass,2);
    
    firstMassIndexInside = find(massesInCollector,1,'first');
    lastMassIndexInside  = find(massesInCollector,1,'last');
    G(iMass, firstMassIndexInside + 1: lastMassIndexInside -1) = deltabeamMassInterp;
    G(iMass, [firstMassIndexInside, lastMassIndexInside]) = deltabeamMassInterp/2;

end

% trim data
hasModelBeam = any(G,2); % magnet masses with beam model mass in collector
G = G(hasModelBeam,:);
data.baseline = data.measPeakIntensity(~hasModelBeam); % baseline when no beam in detector
data.measPeakIntensityBLcorr = data.measPeakIntensity - ...
                               mean(data.baseline(2:end)); % 1st measured intensity is dodgy
data.magnetMasses = data.magnetMasses(hasModelBeam);
data.measPeakIntensity = data.measPeakIntensityBLcorr(hasModelBeam);

% WLS and NNLS
GB = G*splineBasis.B;
Wdata = diag(1./max(data.measPeakIntensity,1));
beamWLS = (GB'*Wdata*GB)\(GB'*Wdata*data.measPeakIntensity);
beamWNNLS = lsqnonneg(chol(Wdata)*GB,chol(Wdata)*data.measPeakIntensity);

% smoothing spline
lambda = 5e-7;
D = diff(eye(beamKnots+bdeg), pord); % 2nd order smoothing, cubic spline;

%Wdata = eye(length(measPeakIntensity));
Gaugmented = [GB; sqrt(lambda)*D];
measAugmented = [data.measPeakIntensity; zeros(beamKnots+bdeg-pord,1)];
wtsAugmented = blkdiag(Wdata, eye(beamKnots+bdeg-pord));
beamPSpline = (Gaugmented'*wtsAugmented*Gaugmented)\(Gaugmented'*wtsAugmented*measAugmented);
[beamNNPspl, resnorm, residual] = lsqnonneg(chol(wtsAugmented)*Gaugmented,chol(wtsAugmented)*measAugmented);





end % function fitPCKData