%% Interpret PeakCentres for beam/peak shape


%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [14, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Mass", "Intensity"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
data = readtable("DVCC18-9 z9 Pb-570-PKC-205Pb-PM-S2B7C1.txt", opts);

% Convert to output type
data = table2array(data);
clear opts


%% plot peak shape

%plot(data(:,1), data(:,2), '-', 'LineWidth', 2)
magnetMassVec = data(:,1);
measIntensityCPS = data(:,2);
measIntensityAnscombe = 2*sqrt(measIntensityCPS+3/8);

%% Input constants

collectorWidthMM = 1;
effectiveRadiusMagnetMM = 540;


%% Calculations

averageMassAMU = mean(magnetMassVec);
collectorWidthAMU = averageMassAMU/effectiveRadiusMagnetMM * collectorWidthMM;

collectorLimits = magnetMassVec + [-collectorWidthAMU, collectorWidthAMU]/2;

nMasses = length(magnetMassVec);

covKernel = zeros(nMasses);
for iMass = 1:nMasses

    covKernel(iMass,:) = (collectorLimits(iMass,1) < magnetMassVec & ...
                          magnetMassVec < collectorLimits(iMass,2) )';

end


%% Solve

% gonna have to use some ridge regression for that kernel

lambda = 0.5;
beam = (covKernel'*covKernel + lambda*eye(nMasses)) \ (covKernel'*measIntensityAnscombe);

% inverse Anscombe transformation
beam = 1/4*beam.^2 - 3/8;
plot(magnetMassVec, beam)