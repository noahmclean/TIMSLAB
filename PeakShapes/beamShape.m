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
Wdata = diag(1./max(measIntensityCPS,1));
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
%plot(magnetMassVec, beam)

% Note: a bad idea because transformation nonlinear?  Probably.  

%% Fix edge effects:
%   1. Extend modeled mass range, 
%   2. Use trapezoidal rule. Better than Simpson's for sharp peaks.

% determine average measured mass spacing. 
% note: last deltaMass is ~10% smaller for unknown reasons, so ignored.
% deltaMagnetMass = mean(magnetMassVec(2:end-1)-magnetMassVec(1:end-2));
% nMassIntervalsInCollector = floor(collectorWidthAMU/deltaMagnetMass);
% nMassIntervalsInCollector = bitor(nMassIntervalsInCollector,1)-1; % even number <= x
% 
% deltaModelMass = collectorWidthAMU/nMassIntervalsInCollector;
% 
% modelMassVec = (magnetMassVec(1):deltaModelMass:magnetMassVec(end))';
% collectorLimits = modelMassVec + [-collectorWidthAMU, collectorWidthAMU]/2;
% modelMassVec = [collectorLimits(1,1); modelMassVec; collectorLimits(end,2)]; 
% nModelMasses = length(modelMassVec);
% nMassIntervals = nModelMasses - 1;
% 
% deltaModelMass = modelMassVec(2:end) - modelMassVec(1:end-1);
% 

%% Starting over: Fix edge effects
%   1. Extend Modeled mass range
%   2. Use Trapezoidal Rule

% determine average measured mass spacing. 
% note: last deltaMass is ~10% smaller for unknown reasons, so ignored.
deltaMagnetMass = mean(magnetMassVec(2:end-1)-magnetMassVec(1:end-2));
nMagnetMasses = length(magnetMassVec);

% mass range that is inside the collector at any measured magnetMass
collectorLimits = magnetMassVec + [-collectorWidthAMU, collectorWidthAMU]/2;

%firstModelMass = magnetMassVec(1); %collectorLimits(1,1);
%lastModelMass = magnetMassVec(end); %collectorLimits(end,2);
firstModelMass = collectorLimits(1,1);
lastModelMass = collectorLimits(end,2);
nModelMasses = ceil((lastModelMass-firstModelMass)/deltaMagnetMass);
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


%% another ridge regression try:


lambda = 1e-12;
beam = (G'*Wdata*G + lambda*eye(nModelMasses)) \ (G'*Wdata*measIntensityCPS);

% inverse Anscombe transformation
%beam = 1/4*beam.^2 - 3/8;
plot(modelMassVec, beam)

% err... pretty bad.

%% Give (better) Tikhonov regularization a try.




%% do some symbolic math
% 
% syms m1 mb m2 y1 yb y2
% yb = y1 + (mb-m1)/(m2-m1) * (y2-y1);
% A1 = (yb+y2)/2 * (m2-mb);
% exp1 = collect(A1,[y1 y2])
% A2 = (y1 + yb)/2 * (mb - m1);
% exp2 = collect(A2, [y1 y2])