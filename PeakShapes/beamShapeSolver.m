%% Interpret PeakCentres for beam/peak shape

% mass spectrometer properties
collectorWidthMM = 1;
effectiveRadiusMagnetMM = 540;
faradayNames = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
ionCounterNames = ["PM", "SEM"];
amplifierResistance = 1e12;

%% Set up the Import Options and import the data

filename = "DVCC18-9 z9 Pb-570-PKC-205Pb-PM-S2B7C1.txt";
%filename = "NBS987 StaticAxH1H2 Bead1Run1-393-PKC-86Sr-Ax-S1B8C1.txt";
data = parsePeakCenterDataFile(filename);
filenameBits = extractBetween(filename, "-", "-");
detectorName = filenameBits(end); clear filenameBits


%% plot peak shape

magnetMasses = data.meas(:,1);
measIntensity = data.meas(:,2);

if any(detectorName == faradayNames) 

    % measure baseline
    
    % subtract baseline

    % use integration time and intensity to come up with some uncertainties


elseif any(detectorName == ionCounterNames)

    % Poisson data, counts = mean
    Wdata = diag(1./max(measIntensity,1));

end


%% collector limits

averageMassAMU = mean(magnetMasses);
collectorWidthAMU = averageMassAMU/effectiveRadiusMagnetMM * collectorWidthMM;

%% setup G

% Starting over: Fix edge effects
%   1. Extend Modeled mass range
%   2. Use Trapezoidal Rule

% determine average measured mass spacing. 
% note: last deltaMass is ~10% smaller for unknown reasons, so ignored.
deltaMagnetMass = mean(magnetMasses(2:end-1)-magnetMasses(1:end-2));
nMagnetMasses = length(magnetMasses);

% mass range that is inside the collector at any measured magnetMass
collectorLimits = magnetMasses + [-collectorWidthAMU, collectorWidthAMU]/2;

%firstModelMass = magnetMassVec(1); %collectorLimits(1,1);
%lastModelMass = magnetMassVec(end); %collectorLimits(end,2);
%firstModelMass = collectorLimits(1,1);
%lastModelMass = collectorLimits(end,2);
firstModelMass = 204.8;
lastModelMass = 205.2;
nModelMasses = ceil(0.95*(lastModelMass-firstModelMass)/deltaMagnetMass);
modelMassVec = linspace(firstModelMass, lastModelMass, nModelMasses);
deltaModelMass = modelMassVec(2)-modelMassVec(1);

% make a new convolution matrix, now with trapezoidal rule and 
G = zeros(nMagnetMasses, nModelMasses);
for iMass = 1:nMagnetMasses % a row for each manget mass

    % massesInCollector are *model* masses
    massesInCollector = collectorLimits(iMass,1) <= modelMassVec & ...
                        modelMassVec <= collectorLimits(iMass,2);
    
    firstMassIndexInside = find(massesInCollector,1,'first');
    lastMassIndexInside  = find(massesInCollector,1,'last');
    G(iMass, firstMassIndexInside + 1: lastMassIndexInside -1) = deltaModelMass;
    G(iMass, [firstMassIndexInside, lastMassIndexInside]) = deltaModelMass/2;

    % add in fractional segments at edges of peak/collector overlap
    if firstMassIndexInside > 1 % if not on front edge
        
        m1 = modelMassVec(firstMassIndexInside-1);
        m2 = modelMassVec(firstMassIndexInside);
        mb = collectorLimits(iMass,1);
        G(iMass, firstMassIndexInside-1) = (m2-mb)^2/(2*(m2-m1));
        G(iMass, firstMassIndexInside) = G(iMass, firstMassIndexInside) + ...
                                    (m2-mb)*(m2-2*m1+mb)/(2*(m2-m1));

    end % if firstMassIndexInside > 1 (not on front edge)

    if lastMassIndexInside < nModelMasses % if not on back edge
        
        m1 = modelMassVec(lastMassIndexInside);
        m2 = modelMassVec(lastMassIndexInside+1);
        mb = collectorLimits(iMass,2);
        G(iMass, lastMassIndexInside) = G(iMass, lastMassIndexInside) + ...
                                (m1 - mb)*(m1-2*m2+mb)/(2*(m2-m1));
        G(iMass, lastMassIndexInside+1) = (mb-m1)^2/(2*(m2-m1));

    end % if firstMassIndexInside > 1 (not on back edge)

end

% note: this is a nicer/more accurate trapezoidal rule G


%% WLS (free and non-negative)

centerBeam1 = (G'*Wdata*G)\(G'*Wdata*measIntensity);
plot(modelMassVec, centerBeam1, '-g', 'LineWidth', 2);

hold on
centerBeam2 = lsqnonneg(chol(Wdata)*G, chol(Wdata)*measIntensity);
plot(modelMassVec, centerBeam2, ':k', 'LineWidth', 2);

dataHat = G*centerBeam2;
figure
hold on
plot(magnetMasses, measIntensity, '-r', 'LineWidth', 2);
plot(magnetMasses, dataHat, ':g', 'LineWidth', 2 )