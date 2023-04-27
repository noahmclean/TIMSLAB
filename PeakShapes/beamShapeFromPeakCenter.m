%% script to fit a beam shape and location to a peak center file
% INPUTS 
% - some information about the mass spectrometer e.g. magnet radius
% - location of the .txt file containing peak center output
% 
% OUTPUTS
% - beam shape as a spline
% - beam width/skewness/kurtosis
% - modeled peak shape to compare with measured peak
%
% for Tripoli project by Noah McLean
% 24 April 2023

%% Setup

addpath("exampleData/")
massSpec = massSpecModel("PhoenixPurdue_ATONA");
%filename = "DVCC18-9 z9 Pb-570-PKC-205Pb-PM-S2B7C1.txt";
%filename = "HY30ZK z10 Pb-1004-PKC-205Pb-PM-S2B7C1.txt";
%filename = "HY30ZK z10 Pb-1004-PKC-207Pb-PM-S4B8C1.TXT";
%filename = "6NHCl dpblank-210204A-169-PKC-208Pb-PM-S5B2C1.txt";
%filename = "NBS987 StaticAxH1H2 Bead1Run1-393-PKC-86Sr-Ax-S1B8C1.txt";
%filename = "Ryan-PKC-270UO-PM-Peak.TXT";
filename = "A519_Pb-2531-PKC-208Pb-H2-S3B6C1.TXT";
data = dataModel(filename);


%% This goes sideways for unclear reasons.

[bestCollectorWidthMM, ~] = fminbnd(@(collectorWidthMM) ...
                                  fitPKCData(massSpec, data, collectorWidthMM), ...
                                  0.5*massSpec.nominalCollectorWidthMM, ...
                                  2.0*massSpec.nominalCollectorWidthMM);

% nWidths = 5000;
% resnormVector = zeros(nWidths,1);
% widthVector = linspace(0.01, 2, nWidths)';
% for iWidth = 1:nWidths
% 
%     resnormi = fitPKCData(massSpec, data, widthVector(iWidth));
%     resnormVector(iWidth) = resnormi;
% 
% end
% plot(widthVector, resnormVector)

[beamNNPspl, resnorm, residual, splineBasis] = fitPKCData(massSpec, data, bestCollectorWidthMM);


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

% with distance as the x-axis
beamMassInterp = splineBasis.beamMassInterp;
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

