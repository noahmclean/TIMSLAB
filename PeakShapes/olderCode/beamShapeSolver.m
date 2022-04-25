%% Interpret PeakCentres for beam/peak shape

% mass spectrometer properties
massSpec = setupMassSpec("PhoenixKansas_1e12");

%% Set up the Import Options and import the data

filename = "DVCC18-9 z9 Pb-570-PKC-205Pb-PM-S2B7C1.txt";
%filename = "NBS987 StaticAxH1H2 Bead1Run1-393-PKC-86Sr-Ax-S1B8C1.txt";
data = parsePeakCenterDataFile(filename);


%% plot peak shape

magnetMasses = data.meas(:,1);
measIntensity = data.meas(:,2);

if any(data.detectorName == massSpec.faradayNames) 

    % measure baseline
    
    % subtract baseline

    % use integration time and intensity to come up with some uncertainties


elseif any(data.detectorName == massSpec.ionCounterNames)

    % Poisson data, counts = mean
    Wdata = diag(1./max(measIntensity,1));

end


%% collector limits

modelMassRange = [204.8 205.2];

[G, modelMasses] = assembleG(magnetMasses, massSpec, modelMassRange);


%% WLS (free and non-negative)

centerBeam1 = (G'*Wdata*G)\(G'*Wdata*measIntensity);
plot(modelMasses, centerBeam1, '-g', 'LineWidth', 2);

hold on
centerBeam2 = lsqnonneg(chol(Wdata)*G, chol(Wdata)*measIntensity);
plot(modelMasses, centerBeam2, ':k', 'LineWidth', 2);

dataHat = G*centerBeam2;
figure
hold on
plot(magnetMasses, measIntensity, '-r', 'LineWidth', 2);
plot(magnetMasses, dataHat, ':g', 'LineWidth', 2 )

