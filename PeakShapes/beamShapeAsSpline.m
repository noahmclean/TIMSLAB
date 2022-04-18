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
collectorWidthAMU = massAtCenterAMU / massSpec.effectiveRadiusMagnetMM * massSpec.collectorWidthMM;
beamWindow = beamWidthAMU*2;
collectorLimits = magnetMasses + [-collectorWidthAMU, collectorWidthAMU]/2;

deltaMagnetMass = magnetMasses(2)-magnetMasses(1);

beamKnots = floor(beamWindow/(deltaMagnetMass)) - 3; % 3 xtra knots for cubic spline
nInterp = 500;

xl = massAtCenterAMU - beamWindow/2;
xr = massAtCenterAMU + beamWindow/2;

beamMassInterp = linspace(xl, xr, nInterp);
B = bbase(beamMassInterp, xl, xr, beamKnots, 3);
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

%% give this a try?
testBeam = [0 0 0 0 0 0 0 0 0 1 4 10 9 8 7 6 5 4 3 2 1 0 0 0 0 0]' * 5*10^5;
GB = G*B;
truePeakIntensity = G*B*testBeam;

figure

subplot(2,2,1)
plot(magnetMasses, truePeakIntensity, 'LineWidth', 2)
subplot(2,2,2)
hold on
plot(beamMassInterp, B*testBeam, 'LineWidth', 2)

beam1 = (G*B)\truePeakIntensity;
plot(beamMassInterp, B*beam1, ':r', 'LineWidth', 2)


% ok!  add noise.

measPeakIntensity = poissrnd(truePeakIntensity*0.2);
Wdata = diag(1./max(measPeakIntensity,1));
beam2 = (GB'*Wdata*GB)\(GB'*Wdata*measPeakIntensity);
beam3 = GB\measPeakIntensity;
subplot(2,2,3)
plot(magnetMasses, measPeakIntensity, 'LineWidth', 2);
subplot(2,2,4)
hold on
plot(beamMassInterp, B*beam3*5, '-b', 'LineWidth', 2)
plot(beamMassInterp, B*beam2*5, ':r', 'LineWidth', 2)

% results not the best, but not horrible!

% now let's try with some smoothing

lambda = 1e-10;
D = diff(eye(beamKnots+3), 2); % 2nd order smoothing, cubic spline;

Gaugmented = [GB; sqrt(lambda)*D];
measAugmented = [measPeakIntensity; zeros(24,1)];
wtsAugmented = blkdiag(Wdata, eye(24));
beamSmooth = (Gaugmented'*wtsAugmented*Gaugmented)\(Gaugmented'*wtsAugmented*measAugmented);
plot(beamMassInterp, B*beamSmooth*5, '-g', 'LineWidth', 2)

% smoother!  need to replace magic numbers (24) above.
