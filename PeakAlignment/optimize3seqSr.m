%% optimize three sequence Sr method as a test/sandbox

%     L3  L2  Ax  H1  H2
% S1          86  87  88
% S2      86  87  88
% S3  86  87  88

% parameterized as 
% m1, m2, m3 axial masses during S1, S2, S3
% d1, d2, d3, d4 distances between L2-L2, L2-Ax, Ax-H1, H1-H2

% constants to use: a,b,c are inverse of dispersion for S1, S2, S3
% note adjusted by ~0.0001-2 for offsets m1, m2, m3
a = 85.9121/540;
b = 86.9031/540;
c = 87.9085/540;

y = repmat([mass.Sr86; mass.Sr87; mass.Sr88], 3, 1);
A = [1 0 0  0  0  0  0;
     1 0 0  0  0  a  0;
     1 0 0  0  0  a  a;
     0 1 0  0 -b  0  0;
     0 1 0  0  0  0  0;
     0 1 0  0  0  b  0;
     0 0 1 -c -c  0  0;
     0 0 1  0 -c  0  0;
     0 0 1  0  0  0  0];
x = A\y;
r = y - A*x;

% results in terms of peak top width
WC = 1.0; % mm, width of collector slit
WB = 0.35; % mm, width of beam with double focusing
Reff = 540; % mm, effective radius of magnet
WT = (WC - WB)*[mass.Sr86 mass.Sr87 mass.Sr88]/540;

isotopeMasses = [mass.Sr86 mass.Sr87 mass.Sr88];
collectorWidths = [0.97 0.97 0.97 0.97 0.97];
nCollectors = length(collectorWidths);

%% reproduce a mass scan using this information

% start simple with a beam shape from beamShape in TIMSLAB
% need: beamDistInterp: x-axis for beam profile, in mm distances
% need: beamShape: fit by beamShape routine

collectorPositions = [-x(4)-x(5) -x(4) 0 x(6) x(6) + x(7)];
collectorEdges = [collectorPositions - 0.5*collectorWidths;
                  collectorPositions + 0.5*collectorWidths];
nCollectorInterpPoints = 3000;
collectorVectors = zeros(nCollectorInterpPoints, nCollectors);
for collectorIdx = 1:5
    collectorVectors(:,collectorIdx) = linspace(collectorEdges(1, collectorIdx), ...
              collectorEdges(2, collectorIdx), nCollectorInterpPoints);
end % for collectorIdx


nScans = 3;
nMassesToScan = 700;

scanLimits = [85.7, 86.12; 86.7, 87.12; 87.7, 88.12];

for iScan = 1:nScans

massScanMass = linspace(scanLimits(iScan, 1), scanLimits(iScan,2), nMassesToScan);
% nMassesToScan = 1;
% massScan = 85.91;

massScanIntensity = zeros(nMassesToScan, nCollectors);



for iMass = 1:nMassesToScan

    axialMass = massScanMass(iMass);
    dispersion = Reff/axialMass;
    isotopePositions = (isotopeMasses - axialMass)*dispersion;

    for collectorIndex = 1:5 % [L3 L2 Ax H1 H2] = [1 2 3 4 5]

        thisCollectorVector = collectorVectors(:, collectorIndex);

        for isotopeIndex = 1:3 % [86Sr 87Sr 88Sr] = [1 2 3]

            beamShapePosition = beamDistInterp + isotopePositions(isotopeIndex);

            beamOverlap = interp1(beamShapePosition, beamShapeNorm, thisCollectorVector, 'pchip', 0);

            massScanIntensity(iMass, collectorIndex) = ...
                massScanIntensity(iMass, collectorIndex) + ...
                trapz(thisCollectorVector, beamOverlap);

        end % for isotopeIndex


    end % for collectorIndex


end % for iMass

subplot(nScans,1,iScan)
plot(massScanMass, massScanIntensity, 'LineWidth', 2)
line([x(iScan) x(iScan)], [0 1.05], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
legend({'L3', 'L2', 'Ax', 'H1', 'H2'})
set(gca, 'FontSize', 16)
ylim([0 1.05])
end % for iScan


