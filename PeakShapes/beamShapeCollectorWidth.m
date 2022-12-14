%% solve for beam shape and collector width using (smoothing) splines

addpath("exampleData/")
massSpec = massSpecModel("PhoenixKansas_1e12");

%filename = "DVCC18-9 z9 Pb-570-PKC-205Pb-PM-S2B7C1.txt";
%filename = "HY30ZK z10 Pb-1004-PKC-205Pb-PM-S2B7C1.txt";
%filename = "HY30ZK z10 Pb-1004-PKC-207Pb-PM-S4B8C1.TXT";
%filename = "6NHCl dpblank-210204A-169-PKC-208Pb-PM-S5B2C1.txt";
filename = "NBS987 StaticAxH1H2 Bead1Run1-393-PKC-86Sr-Ax-S1B8C1.txt";
%filename = "Ryan-PKC-270UO-PM-Peak.TXT";
data = dataModel(filename);

%% 

% calculations about the setup, depend on data and mass spec

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
lambda = 1e-6;
D = diff(eye(beamKnots+bdeg), pord); % 2nd order smoothing, cubic spline;

%Wdata = eye(length(measPeakIntensity));
Gaugmented = [GB; sqrt(lambda)*D];
measAugmented = [data.measPeakIntensity; zeros(beamKnots+bdeg-pord,1)];
wtsAugmented = blkdiag(Wdata, eye(beamKnots+bdeg-pord));
beamPSpline = (Gaugmented'*wtsAugmented*Gaugmented)\(Gaugmented'*wtsAugmented*measAugmented);
beamNNPspl = lsqnonneg(chol(wtsAugmented)*Gaugmented,chol(wtsAugmented)*measAugmented);


%% determine peak width

% note: need a different fit with a negative baseline for Faraday.

beamShape = splineBasis.B*beamNNPspl; % weighted non-negative least squares

[maxBeam, maxBeamIndex] = max(beamShape);
thesholdIntensity = 0.05 * maxBeam;

peakLeft = beamShape(1:maxBeamIndex);
leftAboveTheshold = peakLeft > thesholdIntensity;
leftThesholdChange = leftAboveTheshold(2:end)-leftAboveTheshold(1:end-1);
leftBoundary = find(leftThesholdChange, 1, 'last') + 1;

peakRight = beamShape(maxBeamIndex:end);
rightAboveThreshold = peakRight > thesholdIntensity;
rightThesholdChange = rightAboveThreshold(1:end-1)-rightAboveThreshold(2:end);
rightBoundary = find(rightThesholdChange, 1, 'first') + maxBeamIndex - 1;

measBeamWidthAMU = beamMassInterp(rightBoundary) - beamMassInterp(leftBoundary);
measBeamWidthMM = measBeamWidthAMU * massSpec.effectiveRadiusMagnetMM/data.peakCenterMass;

% with distance as the x-axis
beamDistInterp = beamMassInterp*massSpec.effectiveRadiusMagnetMM/data.peakCenterMass - ...
                 data.peakCenterMass*massSpec.effectiveRadiusMagnetMM/data.peakCenterMass;

% some moments, in distance space (units = mm)
beamShapeNorm = beamShape / trapz(beamDistInterp', beamShape); % normalized 
beamMean = trapz(beamDistInterp', beamDistInterp' .* beamShapeNorm);
beamVariance = trapz(beamDistInterp', (beamDistInterp'-beamMean).^2 .* beamShapeNorm);
beamSkewness = trapz(beamDistInterp', (beamDistInterp'-beamMean).^3 .* beamShapeNorm)...
               / beamVariance^(3/2);
beamKurtosis = trapz(beamDistInterp', (beamDistInterp'-beamMean).^4 .* beamShapeNorm)...
               / beamVariance^(4/2);


%%
figure('Position', [500 500 1200 600]);
subplot(1,2,1); hold on
plot(beamMassInterp, beamShape, '-k', 'LineWidth', 2)
plot(beamMassInterp(leftBoundary), beamShape(leftBoundary),   '.r', 'MarkerSize', 20)
plot(beamMassInterp(rightBoundary), beamShape(rightBoundary), '.r', 'MarkerSize', 20)
line([beamMassInterp(leftBoundary) beamMassInterp(rightBoundary)], ...
     [beamShape(leftBoundary), beamShape(rightBoundary)], ...
     'LineStyle', '--', 'LineWidth', 2, 'Color', 'b')

set(gca, "PlotBoxAspectRatio", [1 1 1])
subplot(1,2,2); hold on
plot(data.magnetMasses, data.measPeakIntensity, '-b', 'LineWidth', 2)
plot(data.magnetMasses, G*beamShape, ':r', 'LineWidth', 2)
set(gca, "PlotBoxAspectRatio", [1 1 1])
