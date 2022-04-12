%% test cases for beam shape model


massSpec = setupMassSpec("PhoenixKansas_1e12");

magnetMasses = 204.65:0.001:205.35;
modelMassRange = [204.75 205.25];

[G, modelMasses] = assembleG(magnetMasses, massSpec, modelMassRange);

nModelMasses = length(modelMasses);
beam = zeros(size(modelMasses));
beamCenter = 

measIntensity = G*beam;
plot(magnetMasses, measIntensity)

xlim([magnetMasses(1) magnetMasses(end)])
