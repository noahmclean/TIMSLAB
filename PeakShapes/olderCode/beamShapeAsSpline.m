%% try a spline basis for this beam shape debacle.

massSpec = setupMassSpec("PhoenixKansas_1e12");

magnetMasses = (204.6:0.01:205.4)';
modelMassRange = [204.85 205.15];
massAtCenterAMU = mean(magnetMasses);
theoreticalBeamWidthMM = 0.35; % mm

%% make a spline

% x = linspace(modelMassRange(1), modelMassRange(2), 500);
% xl = 204.85;
% xr = 205.15;
% nseg = 10;
% bdeg = 3;
% 
% B = bbase(x, xl, xr, nseg, bdeg);
% u = [0 0 1 0 0.5 1.5 1 0.5 0.5 0.2 0.5 0.2 1]';
% y = B*u;

% hold on
% plot(x,y)
% plot(linspace(204.85, 205.15, 13), u)

%% now try again, with limited number of spline segments, no penalty

beamWidthAMU = massAtCenterAMU / massSpec.effectiveRadiusMagnetMM * theoreticalBeamWidthMM;
%collectorWidthAMU = massAtCenterAMU / massSpec.effectiveRadiusMagnetMM * massSpec.collectorWidthMM;
collectorWidthAMU = 0.36;
beamWindow = beamWidthAMU*2;
collectorLimits = magnetMasses + [-collectorWidthAMU, collectorWidthAMU]/2;

deltaMagnetMass = magnetMasses(2)-magnetMasses(1);

bdeg = 3; % order of spline (= order of polynomial pieces)
pord = 2; % order of differences 
beamKnots = 2*floor(beamWindow/(deltaMagnetMass)) - bdeg; % 3 xtra knots for cubic spline
nInterp = 500; % number of interpolated segments

xl = massAtCenterAMU - beamWindow/2;
xr = massAtCenterAMU + beamWindow/2;

beamMassInterp = linspace(xl, xr, nInterp);
B = bbase(beamMassInterp, xl, xr, beamKnots, bdeg);
deltabeamMassInterp = beamMassInterp(2)-beamMassInterp(1);

% discrete convolution setup with design matrix G

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

% %% give this a try?
% testBeam = [0 0 0 0 0 0 0 0 0 1 4 10 9 8 7 6 5 4 3 2 1 0 0 0 0 0]' * 5*10^4;
% GB = G*B;
% truePeakIntensity = G*B*testBeam;
% 
% figure('Position', [1 1 1200 1000], 'Units', 'pixels')
% 
% subplot(2,2,1)
% plot(magnetMasses, truePeakIntensity, 'LineWidth', 2)
% subplot(2,2,2)
% hold on
% plot(beamMassInterp, B*testBeam, 'LineWidth', 2)
% 
% beam1 = (G*B)\truePeakIntensity;
% plot(beamMassInterp, B*beam1, ':r', 'LineWidth', 2)
% 
% 
% % ok!  add noise.
% integrationTime = 0.2; % seconds
% measPeakIntensity = poissrnd(truePeakIntensity*integrationTime);
% Wdata = diag(1./max(measPeakIntensity,1));
% beam2 = (GB'*Wdata*GB)\(GB'*Wdata*measPeakIntensity);
% beam3 = GB\measPeakIntensity;
% subplot(2,2,3)
% plot(magnetMasses, measPeakIntensity, 'LineWidth', 2);
% subplot(2,2,4)
% hold on
% plot(beamMassInterp, B*beam3*5, '-b', 'LineWidth', 2)
% plot(beamMassInterp, B*beam2*5, ':r', 'LineWidth', 2)
% 
% % results not the best, but not horrible!
% 
% % now let's try with some smoothing
% 
% lambda = 1e-10;
% D = diff(eye(beamKnots+bdeg), pord); % 2nd order smoothing, cubic spline;
% 
% Gaugmented = [GB; sqrt(lambda)*D];
% measAugmented = [measPeakIntensity; zeros(beamKnots+bdeg-pord,1)];
% wtsAugmented = blkdiag(Wdata, eye(beamKnots+bdeg-pord));
% beamSmooth = (Gaugmented'*wtsAugmented*Gaugmented)\(Gaugmented'*wtsAugmented*measAugmented);
% plot(beamMassInterp, B*beamSmooth*5, '-g', 'LineWidth', 2)
% 
% % smoother!  need to replace magic numbers (24) above.
% 

%% try on real data

filename = "DVCC18-9 z9 Pb-570-PKC-205Pb-PM-S2B7C1.txt";
%filename = "NBS987 StaticAxH1H2 Bead1Run1-393-PKC-86Sr-Ax-S1B8C1.txt";
data = parsePeakCenterDataFile(filename);

magnetMasses = data.meas(:,1);
measPeakIntensity = data.meas(:,2);
massAtCenterAMU = data.peakCenterMass;

collectorLimits = magnetMasses + [-collectorWidthAMU, collectorWidthAMU]/2;
deltaMagnetMass = magnetMasses(2)-magnetMasses(1);

bdeg = 3; % order of spline (= order of polynomial pieces)
pord = 2; % order of differences 
nInterp = 5000; % number of interpolated segments

xl = massAtCenterAMU - beamWindow/2;
xr = massAtCenterAMU + beamWindow/2;

beamMassInterp = linspace(xl, xr, nInterp);
B = bbase(beamMassInterp, xl, xr, beamKnots, bdeg);
deltabeamMassInterp = beamMassInterp(2)-beamMassInterp(1);

% discrete convolution setup with design matrix G

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

% solution for real data

% trim data
hasModelBeam = any(G,2); % magnet masses with beam model mass in collector
G = G(hasModelBeam,:);
magnetMasses = magnetMasses(hasModelBeam);
measPeakIntensity = measPeakIntensity(hasModelBeam);

GB = G*B;
Wdata = diag(1./max(measPeakIntensity,1));
%Wdata = eye(length(measPeakIntensity));
beam4 = (GB'*Wdata*GB)\(GB'*Wdata*measPeakIntensity);
beam6 = lsqnonneg(chol(Wdata)*GB,chol(Wdata)*measPeakIntensity);


figure
subplot(1,2,1)
plot(beamMassInterp, B*beam4, '-k', 'LineWidth', 2)
subplot(1,2,2); hold on
plot(magnetMasses, measPeakIntensity, '-b', 'LineWidth', 2)
plot(magnetMasses, GB*beam4, ':r', 'LineWidth', 2)

% try some regularization
%%

lambda = 1e-8;
D = diff(eye(beamKnots+bdeg), pord); % 2nd order smoothing, cubic spline;

%Wdata = eye(length(measPeakIntensity));
Gaugmented = [GB; sqrt(lambda)*D];
measAugmented = [measPeakIntensity; zeros(beamKnots+bdeg-pord,1)];
wtsAugmented = blkdiag(Wdata, eye(beamKnots+bdeg-pord));
beam5 = (Gaugmented'*wtsAugmented*Gaugmented)\(Gaugmented'*wtsAugmented*measAugmented);

figure
subplot(1,3,1); hold on
plot(beam5)
plot(beam5, '.', 'MarkerSize', 20)
subplot(1,3,2)
plot(beamMassInterp, B*beam5)
subplot(1,3,3); hold on
plot(magnetMasses, measPeakIntensity, '-k')
plot(magnetMasses, GB*beam5, '-b', 'LineWidth', 2)


%% try to manually improve residuals

GB = G*B;
Wdata = diag(1./max(measPeakIntensity,1));
%Wdata = eye(length(measPeakIntensity));
beam4 = (GB'*Wdata*GB)\(GB'*Wdata*measPeakIntensity);
% figure
% subplot(1,2,1)
% plot(beamMassInterp, B*beam4, '-k', 'LineWidth', 2)
% subplot(1,2,2); hold on
% plot(magnetMasses, measPeakIntensity, '-b', 'LineWidth', 2)
% plot(magnetMasses, GB*beam4, ':r', 'LineWidth', 2)

r = GB*beam4 - measPeakIntensity;
r2 = (GB*beam4 - measPeakIntensity).^2;
chi2 = r2 .* diag(Wdata);

