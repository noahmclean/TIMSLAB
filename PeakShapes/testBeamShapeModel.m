%% test cases for beam shape model


massSpec = setupMassSpec("PhoenixKansas_1e12");

magnetMasses = 204.6:0.01:205.4;
modelMassRange = [204.85 205.15];

[G, modelMasses] = assembleG(magnetMasses, massSpec, modelMassRange);

nModelMasses = length(modelMasses);
deltaModelMass = modelMasses(2) - modelMasses(1);

beam = zeros(size(modelMasses));
if ~mod(nModelMasses,2) % if even
    beamCenterIndex = [nModelMasses/2 nModelMasses/2+1];
else
    beamCenterIndex = ceil(nModelMasses/2);
end % if

beamWidthAMU = 205/massSpec.effectiveRadiusMagnetMM*0.30;
beamDensity = 1e5; 
halfBeamWidthIndices = round(beamWidthAMU/deltaModelMass/2);

beamStartIndex = beamCenterIndex(1) - halfBeamWidthIndices;
beamStopIndex  = beamCenterIndex(end) + halfBeamWidthIndices;
beam(beamStartIndex:beamStopIndex) = beamDensity;

trueCPS = G*beam;
trueCounts = G*beam * 0.2; % 0.2 seconds of integration
measCounts = poissrnd(trueCounts);

plot(magnetMasses, trueCounts, magnetMasses, measCounts)
xlim([magnetMasses(1) magnetMasses(end)])

Wdata = diag(1./max(measCounts,1));
beamSolution0 = G\trueCounts;
beamSolution1 = (G'*Wdata*G)\(G'*Wdata*measCounts);
beamSolution2 = lsqnonneg(chol(Wdata)*G,chol(Wdata)*measCounts);

meas1 = G*beamSolution1;
meas2 = G*beamSolution2;

subplot(2,1,1)
plot(modelMasses, beamSolution0, modelMasses, beamSolution1)
subplot(2,1,2)
plot(modelMasses, beamSolution0, modelMasses, beamSolution2)


%% test total beam

% this model assumes a piecewise linear beam shape, so if beam
% is a square, the trapezoidal rule as formulated will assume a 
% triangular beam segment at each edge, where beam ramps from
% zero to beamDensity.  The sum of the two triangles makes for 
% one extra interval of deltaModelMass

% total beams should be same if G is correctly integrating
totalBeam = (modelMasses(beamStopIndex) - ...
             modelMasses(beamStartIndex) + 1)*beamDensity;
maxMeasuredBeamIntensity = max(trueCPS);

