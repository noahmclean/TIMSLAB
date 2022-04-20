%% solve for beam shape and collector width using (smoothing) splines

addpath("exampleData/")
massSpec = setupMassSpec("PhoenixKansas_1e12");
%filename = "DVCC18-9 z9 Pb-570-PKC-205Pb-PM-S2B7C1.txt";
%filename = "HY30ZK z10 Pb-1004-PKC-205Pb-PM-S2B7C1.txt";
filename = "HY30ZK z10 Pb-1004-PKC-207Pb-PM-S4B8C1.TXT";
%filename = "6NHCl dpblank-210204A-169-PKC-208Pb-PM-S5B2C1.txt";
%filename = "NBS987 StaticAxH1H2 Bead1Run1-393-PKC-86Sr-Ax-S1B8C1.txt";

data = parsePeakCenterDataFile(filename);

%%

theoreticalBeamWidthMM = 0.35; % mm
%collectorWidths = 0.951422845691383; %linspace(0.9, 1.1, 500);

collectorWidths = 0.95135;

chi2 = zeros(size(collectorWidths));
for i = 1:length(collectorWidths)

collectorWidthMM = collectorWidths(i);

massAtCenterAMU = data.peakCenterMass;
collectorWidthAMU = massAtCenterAMU / massSpec.effectiveRadiusMagnetMM * collectorWidthMM;
beamWidthAMU      = massAtCenterAMU / massSpec.effectiveRadiusMagnetMM * theoreticalBeamWidthMM;

magnetMasses = data.meas(:,1);
measPeakIntensity = data.meas(:,2);

collectorLimits = magnetMasses + [-collectorWidthAMU, collectorWidthAMU]/2;
deltaMagnetMass = magnetMasses(2)-magnetMasses(1);

beamWindow = beamWidthAMU*2;

bdeg = 3; % order of spline (= order of polynomial pieces)
pord = 2; % order of differences 
beamKnots = ceil(beamWindow/(deltaMagnetMass)) - 2*bdeg; % 3 xtra knots for cubic spline
nInterp = 1000; % number of interpolated segments

xl = massAtCenterAMU - beamWindow/2;
xr = massAtCenterAMU + beamWindow/2;

beamMassInterp = linspace(xl, xr, nInterp);
B = bbase(beamMassInterp, xl, xr, beamKnots, bdeg);
deltabeamMassInterp = beamMassInterp(2)-beamMassInterp(1);

nMagnetMasses = length(magnetMasses);
G = zeros(nMagnetMasses, nInterp);
for iMass = 1:nMagnetMasses % a row for each manget mass

    % massesInCollector are *model* masses
    massesInCollector = collectorLimits(iMass,1) <= beamMassInterp & ...
                        beamMassInterp <= collectorLimits(iMass,2);
    
    firstMassIndexInside = find(massesInCollector,1,'first');
    lastMassIndexInside  = find(massesInCollector,1,'last');
    G(iMass, firstMassIndexInside + 1: lastMassIndexInside -1) = deltabeamMassInterp;
    G(iMass, [firstMassIndexInside, lastMassIndexInside]) = deltabeamMassInterp/2;

end

% trim data
hasModelBeam = any(G,2); % magnet masses with beam model mass in collector
G = G(hasModelBeam,:);
magnetMasses = magnetMasses(hasModelBeam);
measPeakIntensity = measPeakIntensity(hasModelBeam);

% WLS and NNLS
GB = G*B;
Wdata = diag(1./max(measPeakIntensity,1));
beamWLS = (GB'*Wdata*GB)\(GB'*Wdata*measPeakIntensity);
beamWNNLS = lsqnonneg(chol(Wdata)*GB,chol(Wdata)*measPeakIntensity);

% smoothing spline
lambda = 1e-6;
D = diff(eye(beamKnots+bdeg), pord); % 2nd order smoothing, cubic spline;

%Wdata = eye(length(measPeakIntensity));
Gaugmented = [GB; sqrt(lambda)*D];
measAugmented = [measPeakIntensity; zeros(beamKnots+bdeg-pord,1)];
wtsAugmented = blkdiag(Wdata, eye(beamKnots+bdeg-pord));
beamPSpline = (Gaugmented'*wtsAugmented*Gaugmented)\(Gaugmented'*wtsAugmented*measAugmented);
beamNNPSspl = lsqnonneg(chol(wtsAugmented)*Gaugmented,chol(wtsAugmented)*measAugmented);

chi2(i) = sum( (GB*beamWNNLS - measPeakIntensity).^2 .* diag(Wdata) );

end % for

[minchi2, minchi2indx] = min(chi2);
bestCollectorWidth = collectorWidths(minchi2indx);

%% determine peak width

beamShape = B*beamNNPSspl;
[maxBeam, maxBeamIndex] = max(beamShape);
thesholdIntensity = 0.02 * maxBeam;

peakLeft = beamShape(1:maxBeamIndex);
leftAboveTheshold = peakLeft > thesholdIntensity;
leftThesholdChange = leftAboveTheshold(2:end)-leftAboveTheshold(1:end-1);
leftBoundary = find(leftThesholdChange, 1, 'last') + 1;

peakRight = beamShape(maxBeamIndex:end);
rightAboveThreshold = peakRight > thesholdIntensity;
rightThesholdChange = rightAboveThreshold(1:end-1)-rightAboveThreshold(2:end);
rightBoundary = find(rightThesholdChange, 1, 'first') + maxBeamIndex - 1;

measBeamWidthAMU = beamMassInterp(rightBoundary) - beamMassInterp(leftBoundary);
measBeamWidthMM = measBeamWidthAMU * massSpec.effectiveRadiusMagnetMM/massAtCenterAMU;


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
plot(magnetMasses, measPeakIntensity, '-b', 'LineWidth', 2)
plot(magnetMasses, GB*beamWNNLS, ':r', 'LineWidth', 2)
set(gca, "PlotBoxAspectRatio", [1 1 1])
