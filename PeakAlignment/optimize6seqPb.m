%% optimize three sequence Sr method as a test/sandbox

%     Ax   H1   H2   H3   H4
% S1            204       206
% S2       204       206  207
% S3  204       206  207  208
% S4       206  207  208
% S5  206  207  208
% S6  207  208

% parameterized as 
% m1, m2, m3 m4 m5 m6: axial masses during S1, S2, S3
% d1, d2, d3, d4 distances between Ax-H1, H1-H2, H2-H3, H3-H4

% constants to use: a,b,c,d,e,f are inverse of dispersion for S1-S6
% note adjusted by ~0.0001-2 for offsets m1, m2, m3
a = 202.0065/540;
b = 202.9909/540;
c = 203.9780/540;
d = 204.9739/540;
e = 205.9722/540;
f = 206.9733/540;

% y is true masses of the isotopes
y = [mass.Pb204; mass.Pb206; mass.Pb204; mass.Pb206; mass.Pb207; % S1-S2
     mass.Pb204; mass.Pb206; mass.Pb207; mass.Pb208; % S3
     mass.Pb206; mass.Pb207; mass.Pb208; % S4
     mass.Pb206; mass.Pb207; mass.Pb208; mass.Pb207; mass.Pb208]; % S5-S6

mBlock = [1 0 0 0 0 0; % S1: m1
          1 0 0 0 0 0; % S1: m1
          0 1 0 0 0 0; % S2: m2
          0 1 0 0 0 0; % S2: m2
          0 1 0 0 0 0; % S2: m2
          0 0 1 0 0 0; % S3: m3
          0 0 1 0 0 0; % S3: m3
          0 0 1 0 0 0; % S3: m3
          0 0 1 0 0 0; % S3: m3
          0 0 0 1 0 0; % S4: m4
          0 0 0 1 0 0; % S4: m4
          0 0 0 1 0 0; % S4: m4
          0 0 0 0 1 0; % S5: m5
          0 0 0 0 1 0; % S5: m5
          0 0 0 0 1 0; % S5: m5
          0 0 0 0 0 1; % S6: m6
          0 0 0 0 0 1];% S6: m6

dBlock = [a a 0 0; % S1: 204
          a a a a; % S1: 206
          b 0 0 0; % S2: 204
          b b b 0; % S2: 206
          b b b b; % S2: 207
          0 0 0 0; % S3: 204
          c c 0 0; % S3: 206
          c c c 0; % S3: 207
          c c c c; % S3: 208
          d 0 0 0; % S4: 206
          d d 0 0; % S4: 207
          d d d 0; % S4: 208
          0 0 0 0; % S5: 206
          e 0 0 0; % S5: 207
          e e 0 0; % S5: 208
          0 0 0 0; % S6: 207
          f 0 0 0];% S6: 208
A = [mBlock dBlock];

x = A\y;
r = y - A*x;
yhat = A*x;

% results in terms of peak top width
WC = 1.0; % mm, width of collector slit
WB = 0.35; % mm, width of beam with double focusing
Reff = 540; % mm, effective radius of magnet
WT = (WC - WB)*[mass.Pb204 mass.Pb206 mass.Pb207 mass.Pb208]/540;

isotopeMasses = [mass.Pb204 mass.Pb206 mass.Pb207 mass.Pb208];
collectorWidths = [0.97 0.97 0.97 0.97 0.97];
nCollectors = length(collectorWidths);

%% reproduce a mass scan using this information

% start simple with a beam shape from beamShape in TIMSLAB
% required variables from beamShape code:
% need: beamDistInterp: x-axis for beam profile, in mm distances
% need: beamShape: fit by beamShape routine

collectorPositions = [0 x(7) sum(x(7:8)) sum(x(7:9)) sum(x(7:10))];
collectorEdges = [collectorPositions - 0.5*collectorWidths;
                  collectorPositions + 0.5*collectorWidths];
nCollectorInterpPoints = 3000;
collectorVectors = zeros(nCollectorInterpPoints, nCollectors);
for collectorIdx = 1:5
    collectorVectors(:,collectorIdx) = linspace(collectorEdges(1, collectorIdx), ...
              collectorEdges(2, collectorIdx), nCollectorInterpPoints);
end % for collectorIdx


% cmap = [50,  205, 50 ;      % lime green, L3
%         210, 105, 30 ;      % chocolate, L2
%         139, 0,   0  ;      % dark red, Ax
%         85,  107, 47 ;      % dark olive green, H1
%         139, 0,   139]/255; % dark magenta, H2

cmap = lines(5);

colororder(cmap);

nScans = 6;
nMassesToScan = 700;

scanLimits = [201.7 202.3;
              202.7 203.3;
              203.7 204.3
              204.7 205.3
              205.7 206.3
              206.7 207.3];

for iScan = 1:nScans

massScanMass = linspace(scanLimits(iScan, 1), scanLimits(iScan,2), nMassesToScan);
% nMassesToScan = 1;
% massScan = 85.91;

massScanIntensity = zeros(nMassesToScan, nCollectors);



for iMass = 1:nMassesToScan

    axialMass = massScanMass(iMass);
    dispersion = Reff/axialMass;
    isotopePositions = (isotopeMasses - axialMass)*dispersion;

    for collectorIndex = 1:5 % [Ax H1 H2 H3 H4] = [1 2 3 4 5]

        thisCollectorVector = collectorVectors(:, collectorIndex);

        for isotopeIndex = 1:4 % [204Pb 206Pb 207Pb 208Pb] = [1 2 3 4]

            beamShapePosition = beamDistInterp + isotopePositions(isotopeIndex);

            beamOverlap = interp1(beamShapePosition, beamShapeNorm, thisCollectorVector, 'pchip', 0);

            massScanIntensity(iMass, collectorIndex) = ...
                massScanIntensity(iMass, collectorIndex) + ...
                trapz(thisCollectorVector, beamOverlap);

        end % for isotopeIndex


    end % for collectorIndex


end % for iMass

%subplot(nScans,1,iScan); hold on
subplot(2,3,iScan); hold on
plot(massScanMass, massScanIntensity, 'LineWidth', 2)
line([x(iScan) x(iScan)], [0 1.05], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')

% for iPeak = 1:3
%     peakCenterMass = yhat(3*iScan-3 + iPeak);
%     peakCenterColor = cmap(-iScan+3 + iPeak,:);
%     plot(peakCenterMass, 1.1, '|', 'MarkerSize', 30, ...
%         'MarkerEdgeColor', peakCenterColor);
% end

legend({'Ax', 'H1', 'H2', 'H3', 'H4'})
set(gca, 'FontSize', 16)
ylim([0 1.1])
end % for iScan

%% Create peak center 'PKC' table
% masses when you peak center isotopes in cups

seqTableMass = ...
    [0          0          mass.Pb204 0          mass.Pb206;
     0          mass.Pb204 0          mass.Pb206 mass.Pb207;
     mass.Pb204 0          mass.Pb206 mass.Pb207 mass.Pb208;
     0          mass.Pb206 mass.Pb207 mass.Pb208 0;
     mass.Pb206 mass.Pb207 mass.Pb208 0          0;
     mass.Pb207 mass.Pb208 0          0          0];

PKC = zeros(size(seqTableMass));
offset = zeros(size(seqTableMass));
offsetIndex = 0;
for iSeq = 1:6
    for jColl = 1:5

        if seqTableMass(iSeq,jColl) > 0 % if peak present
        PKC(iSeq,jColl) = seqTableMass(iSeq,jColl) / ...
                          (1 + collectorPositions(jColl)/540);
        
        offsetIndex = offsetIndex + 1;
        offset(iSeq,jColl) = -r(offsetIndex);

        end % if peak present

    end % for jColl
end % for iSeq


