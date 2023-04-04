%% make synthetic .txt data file for a Phoenix TIMS
% use real method files parsed by FaradayRelativeEfficiency codes
% use a mass spectrometer model from PeakShapes class definition
% add an ion beam intensity time-series
% output: .TXT file in syndata subfolder
% 
% created by Noah McLean for Tripoli and FaradayRelativeEfficiencies 
% on 28-Feb-2023

%% 1. Setup 

% % create a new sample (or use a reference material)
% s.name    = "NBS981Measurement";
% s.element = "Pb";
% s.species =            ["204Pb", "205Pb", "206Pb", "207Pb", "208Pb"];
% s.relativeAbundances = [0.0590074, 1e-6,   1,       0.914683, 2.1681];
% spl = sample(s.name, s.element, s.species, s.relativeAbundances);
% methodName = "Pb 4-5-6-7-8 Daly 10-5-5-5-2 sec.TIMSAM";

s.name    = "PbTwoIsotopeTwoSequence_1Mcps";
s.element = "Pb";
s.species =            ["206Pb", "208Pb"];
s.relativeAbundances = [1, 1+eps];
spl = sample(s.name, s.element, s.species, s.relativeAbundances);
methodName = "Pb TwoIsotopeTwoSeq 206-8 Ax-PM-H1.TIMSAM";

% name the data file -- refactor?
synDataFileName = s.name;

% add TIMSLAB to path to use its functions/classes
addpath(genpath("../../../TIMSLAB"));
massSpec = massSpecModel("PhoenixKansas_1e12");
nBlocks = 10;

intensityFunction = @(t) 1e6*ones(size(t)); % cps of major isotope
%intensityFunction = @(ampl, freq, minInt, t)  ...
%      ampl*(ceil(freq*t)-freq*t)+minInt; % sawtooth

% isotopic fractionation for Faradays and Ion Counters
betaFaraday = @(t) -0.20*ones(size(t)); % 0.10%/amu at Pb mass
betaDaly    = @(t) -0.32*ones(size(t)); % 0.15%/amu at Pb mass
% using (a/b)meas = (a/b)true*(Ma/Mb)^beta

% collector relative efficiencies
%          PM  RS L5  L4  L3  L2  Ax  H1  H2  H3  H4
CREtrue = @(t) [0.9 1  1   1   1   1   1   1   1   1   1].*ones(size(t));

% baselines
darkNoise = [0.2 0]; % cps, [PM RS]
%refVolts:   L5    L4    L3    L2    Ax     H1    H2    H3    H4
refVolts  = @(t) [-1e-2 -2e-2 -1e-2 -2e-2 -3e-2 -1e-2 -2e-2 -1e-2 -2e-2].*ones(size(t));
kB = 1.38064852e-23;
tempInK = 290;

method = parseTIMSAM(methodName);
% CollNames matches headers in data file for BL and OP
CollNames = ["PM", "RS", "L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
method = processMethod(method, CollNames);

% timing
tStart = 150; % seconds, time of first baseline measurement
tBetweenBlocks = 200;


%% 2. Piece out timing from method

[nBaselines, nOnPeaks, nCyclesPerBlock, integrations, settleTime] ...
                                                               = getMethodTiming(method);

%% 3. Write the header

header = ...
    ["#HEADER";
     "Analysis";
     "Version,"        + "2.0.11,1.04" % consistent with v.2.0.11 of Isolynx, 2022
     "Filename,"       + synDataFileName;
     "MethodName,"     + methodName;
     "MethodPath,"     + "/TIMSLAB/SyntheticData/methodFiles";
     "IsoWorksMethod," + "This is synetic data from the TIMSLAB repository";
     "FolderPath,"     + "/TIMSLAB/SyntheticData/syndata";
     "Corrected,"      + "Yes"; % Presumably corrected for collector gain and efficiency
     "BChannels,"      + "No"; % no ATONA BChannels yet
     "TimeZero,"       + string(datetime("now", "format", "d MMMM yyyy HH:mm:ss.SSS"));
     "";
     "#COLLECTORS";
     "Name,Type,Resistor,Gain,Efficiency,DT"];

writematrix(header, "../syndata/"+synDataFileName, "QuoteStrings", "none")


%% 4. Write the collector block

%fid = fopen("../syndata/" + synDataFileName, 'wt+');
nIonCounters = size(massSpec.ionCounterNames,2);
nFaradays = size(massSpec.faradayNames,2);
nCollectors = nIonCounters + nFaradays;
collBlock = strings(nCollectors, 6);

% ion counters
collBlock(1:nIonCounters,1)     = massSpec.ionCounterNames';
collBlock(1:nIonCounters,2)     = massSpec.ionCounterTypes';
collBlock(1:nIonCounters,3)     = compose("%1.0e", 1e11);
collBlock(1:nIonCounters,4:5)   = compose("%1.9f", 1); % gains & efficiences = 1
collBlock(1:nIonCounters,6)     = compose("%1.4f", massSpec.ionCounterDeadTimes');

% faradays
fStart = nIonCounters + 1; % start index for faradays
collBlock(fStart:end,1)   = massSpec.faradayNames';
collBlock(fStart:end,2)   = "F";
collBlock(fStart:end,3)   = compose("%1.0e", massSpec.amplifierResistance');
collBlock(fStart:end,4:5) = compose("%1.9f", 1);
collBlock(fStart:end,6)   = compose("%1.4f", 0);

collBlock = strtrim(collBlock); % remove leading and trailing whitespace
writematrix(collBlock, "../syndata/"+synDataFileName, ...
    "QuoteStrings", "none", "WriteMode", "append")

% blank line the lazy way
writematrix("", "../syndata/"+synDataFileName, ...
    "QuoteStrings", "none", "WriteMode", "append")


%% 5. Simulate block data, save off in matrices

% make a string array for BL and OP data (serial columns and intensities)
totalIntegrationsBL = sum(integrations.BL.n)*nCyclesPerBlock*nBlocks;
totalIntegrationsOP = sum(integrations.OP.n)*nCyclesPerBlock*nBlocks;
nSerialColumns = 7;
nCollectors = nIonCounters + nFaradays;
ionCounterColumnIndices = nSerialColumns+1:nSerialColumns+nIonCounters;
faradayColumnIndices = nSerialColumns+nIonCounters+1:nSerialColumns+nCollectors;
ionCounterMethodIndices = 1:nIonCounters; % indices of ion counters in method/F_ind
faradayMethodIndices = nIonCounters + 1: nCollectors; % indices of faradays in F_ind

logMassRatios = spl.logNormMasses;
logAbunRatios = spl.logRatioAbundances;
massVector = spl.massVector;
denominatorMass = massVector(spl.denominatorIsotopeIndex);


BL = [];
OP = [];

tCurrent = tStart;
for iBlock = 1:nBlocks

    % baselines first
    for iBLseq = 1:length(integrations.BL.n)

        % make a string array BLseq to contain data from this BL sequence
        nIntegrations = integrations.BL.n(iBLseq);
        integrationPeriod = integrations.BL.integPeriods(iBLseq);
        BLseq = strings(nIntegrations, nSerialColumns + nCollectors);
        
        % populate serial columns
        BLseq(:,1) = "BL" + iBLseq; % ID
        BLseq(:,2) = iBlock; % Block
        BLseq(:,3) = 0; % Cycle
        BLseq(:,4) = (1:nIntegrations)'; % Integ
        BLseq(:,5) = "BL" + iBLseq; % PeakID
        BLseq(:,6) =  num2str(method.axialMasses.BL(iBLseq), '%1.4f');
        
        % times for each integration; update tCurrent
        tCurrent = tCurrent + settleTime.BL(iBLseq);
        tStop = tCurrent + (nIntegrations-1)*integrationPeriod;
        tvector = linspace(tCurrent, tStop, nIntegrations)';
        BLseq(:,7) = num2str(tvector, '%1.7f');
        tCurrent = tStop + integrationPeriod;

        % ion counter baselines (dark noise), units of cps
        lambda = repmat(darkNoise*integrationPeriod, nIntegrations,1);
        darkNoiseIC = random('poisson', lambda);
        BLseq(:,ionCounterColumnIndices) = compose("%1.12e", darkNoiseIC);

        % faraday baselines (Johnson noise), units of volts
        s2 = 4 * kB * tempInK * massSpec.amplifierResistance / integrationPeriod;
        mu = refVolts(tvector);
        sigma = repmat(sqrt(s2), nIntegrations, 1);
        johnsonNoiseF = random('normal', mu, sigma);
        BLseq(:,faradayColumnIndices) = compose("%1.12e", johnsonNoiseF);

        % update BL with this sequence
        BL = [BL; BLseq]; %#ok<AGROW>

    end % for iBLseq
    
    for iCycle = 1:nCyclesPerBlock
    % onpeaks next
    for iOPseq = 1:length(integrations.OP.n)

        % make a string array OPseq to contain data from this OP sequence
        nIntegrations = integrations.OP.n(iOPseq);
        integrationPeriod = integrations.OP.integPeriods(iOPseq);
        OPseq = strings(nIntegrations, nSerialColumns + nCollectors);

        % populate serial columns
        OPseq(:,1) = "OP" + iOPseq; % ID
        OPseq(:,2) = iBlock; % Block
        OPseq(:,3) = iCycle; % Cycle
        OPseq(:,4) = (1:nIntegrations)'; % Integ
        OPseq(:,5) = string(method.onpeaks(iOPseq).MassID); % PeakID
        OPseq(:,6) =  num2str(method.axialMasses.OP(iOPseq), '%1.4f');

        % times for each integration; update tCurrent
        tCurrent = tCurrent + settleTime.OP(iOPseq); % settle before 1st integration
        tStop = tCurrent + (nIntegrations-1)*integrationPeriod;
        tvector = linspace(tCurrent, tStop, nIntegrations)';
        OPseq(:,7) = num2str(tvector, '%1.7f');
        tCurrent = tStop + integrationPeriod;

        % on peaks:
        % 1a. dark noise
        lambda = repmat(darkNoise*integrationPeriod, nIntegrations,1);
        darkNoiseMatrix = random('poisson', lambda);
        % 1b. Johnson noise
        s2 = 4 * kB * tempInK * massSpec.amplifierResistance / integrationPeriod;
        mu = refVolts(tvector);
        sigma = repmat(sqrt(s2), nIntegrations, 1);
        johnsonNoiseMatrix = random('normal', mu, sigma);
        % 1c. Concatenate zero-intensity noise sources
        noiseMatrixSeq = [darkNoiseMatrix johnsonNoiseMatrix];
        
        % 2. distribute isotope parameters into matrices with columns defined by collector array
        % syntax is collectorMatrix(collectorRefsInMethod) = isotopeMatrix(collectorIndicesUsed)
        % simplified, illustrative example of syntax above:
        % x = [1 3 5 7 9]                 vector of scaled intensities
        % y = [0 0 0 0 2 0 0 3 0 0 4 0 5] row of F_ind
        % z = zeros(1,13)                 preallocate data matrix
        % z(find(y>0)) = x(y(y>0))        distribute elements of x to positions in y
        collectorRefsForSequence = method.F_ind(iOPseq,:);
        collectorRefsInMethod = find(collectorRefsForSequence > 0);
        collectorIndicesUsed = collectorRefsForSequence(collectorRefsForSequence>0);

        % 3. betas for ion counters and faradays
        betaMatrixSeq = ones(nIntegrations,nCollectors); % 1 * -Inf = -Inf, exp(-Inf) = 0
        betaMatrixDaly    = repmat(betaDaly(tvector),    1, nIonCounters);
        betaMatrixFaraday = repmat(betaFaraday(tvector), 1, nFaradays);
        betaMatrixAll = [betaMatrixDaly betaMatrixFaraday];
        betaMatrixSeq(:,collectorRefsInMethod) = betaMatrixAll(:,collectorIndicesUsed);
        
        % ion beam intensities (true, logged true)
        speciesIntensities = intensityFunction(tvector)*spl.relativeAbundances;
        trueIntSeq = zeros(nIntegrations, nCollectors);
        trueIntSeq(:,collectorRefsInMethod) = speciesIntensities(:,collectorIndicesUsed);
        logTrueIntSeq = log(trueIntSeq); 
        
        % isotope mass ratios (for fractionation correction)
        logIsotopeMasses = repmat(spl.logNormMasses, nIntegrations, 1);
        logIsoMassSeq = -inf(nIntegrations, nCollectors);
        logIsoMassSeq(:,collectorRefsInMethod) = logIsotopeMasses(:,collectorIndicesUsed);
        
        % collector efficiencies
        CREtrueSeq = CREtrue(tvector);
        collectEffSeq = -inf(nIntegrations, nCollectors);
        collectEffSeq(:,collectorRefsInMethod) = CREtrueSeq(:,collectorIndicesUsed);

        % no tails yet
        tailsSeq = zeros(nIntegrations, nCollectors);

        % 4. put it all together, one step at a time
        % 4a. intensities, isotopically fractionated
        intFractSeq = exp(logTrueIntSeq + betaMatrixSeq .* logIsoMassSeq) + tailsSeq;

        % 4b. apply gain to ion counters only (for use in shot noise calculation)
        intICEffFractSeq = intFractSeq;
        intICEffFractSeq(:,ionCounterMethodIndices) = ...
            intFractSeq(:,ionCounterMethodIndices) .* CREtrueSeq(:,ionCounterMethodIndices);
        
        % 4c. undo dead time correction to get measured counts on ion counter
        dtUnCorrFactor = 1 ./ (1 + repmat(massSpec.ionCounterDeadTimes*1e-9, nIntegrations,1) .* ...
              intICEffFractSeq(:,ionCounterMethodIndices) ); % un-correction factor for dt
        intIonArrivals = intICEffFractSeq;
        intIonArrivals(:,ionCounterMethodIndices) = ...
            intFractSeq(:,ionCounterMethodIndices) .* dtUnCorrFactor;

        % 4d. Add up true intensity/voltage contributors. Mass spec model:
        % i_a = ( exp( log(a/b) + log(i_b) + beta*log(Ma/Mb) ) + tail_a*i_b ) ...
        %       * exp(log(ea/eAx) + refVolts 
        % for i_a is intensity of isotope a, Ma mass of a, ea efficiency
        % for cup with a, tail_a is peak tail contribution to a, beta is fractionation
        
        intPlusShotNoise = random('poisson', intIonArrivals);

        % 4d. redo dead time correction, as implemented in Isolynx software
        dtCorrFactor = 1 ./ (1 - repmat(massSpec.ionCounterDeadTimes*1e-9, nIntegrations,1) .* ...
              intPlusShotNoise(:,ionCounterMethodIndices) ); % correction factor for dt
        intPlusShotNoise(:,ionCounterMethodIndices) = ...
            intPlusShotNoise(:,ionCounterMethodIndices) .* dtCorrFactor;

        % 4e. apply gain to faradays as well, convert to volts
        intPlusShotNoise(:,faradayMethodIndices) = ...
        intPlusShotNoise(:,faradayMethodIndices) .* CREtrueSeq(:,faradayMethodIndices);
        intPlusShotNoise(:,faradayMethodIndices) = ...
        intPlusShotNoise(:,faradayMethodIndices) .* repmat(massSpec.voltsPerCPS, nIntegrations, 1);

        % 4f. add in Johnson and dark noise, plus refVolts for Faradays
        outputOP = intPlusShotNoise + noiseMatrixSeq;
        OPseq(:,[ionCounterColumnIndices faradayColumnIndices]) = compose("%1.12e", outputOP);

        % update OP with this sequence
        OP = [OP; OPseq]; %#ok<AGROW>

    end % for iOPseq

    tCurrent = tCurrent + settleTime.flyBack;

    end % for iCycle

    tCurrent = tCurrent + tBetweenBlocks;


end % for iBlock = 1:nBlocks


%% 6. Write baselines to file

writematrix("#BASELINES", "../syndata/"+synDataFileName, ...
    "QuoteStrings", "none", "WriteMode", "append")

writematrix(["ID","Block","Cycle","Integ","PeakID","AxMass","Time",  ...
            massSpec.ionCounterNames, ...
            massSpec.faradayNames], "../syndata/"+synDataFileName, ...
            "QuoteStrings", "none", "WriteMode", "append")

% add some spaces after commas in BL
spaceArray = strings(size(BL));
spaceArray(:,2:end) = " ";
BL = spaceArray + BL;

writematrix(BL, "../syndata/"+synDataFileName, ...
    "QuoteStrings", "none", "WriteMode", "append")

writematrix(" ", "../syndata/"+synDataFileName, ...
    "QuoteStrings", "none", "WriteMode", "append")


%% 7. Write onPeaks to file

writematrix("#ONPEAK", "../syndata/"+synDataFileName, ...
    "QuoteStrings", "none", "WriteMode", "append")

writematrix(["ID","Block","Cycle","Integ","PeakID","AxMass","Time",  ...
            massSpec.ionCounterNames, ...
            massSpec.faradayNames], "../syndata/"+synDataFileName, ...
            "QuoteStrings", "none", "WriteMode", "append")

% add some spaces after commas in BL
spaceArray = strings(size(OP));
spaceArray(:,2:end) = " ";
OP = spaceArray + OP;

writematrix(OP, "../syndata/"+synDataFileName, ...
    "QuoteStrings", "none", "WriteMode", "append")

writematrix([" "; "#END,AnalysisCompleted"], "../syndata/"+synDataFileName, ...
    "QuoteStrings", "none", "WriteMode", "append")



%% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%

function [nBaselines, nOnPeaks, nCyclesPerBlock, integrations, settleTime] = ...
                                                              getMethodTiming(method)
nBaselines = size(method.baselines,2);
nOnPeaks   = size(method.onpeaks,2);

% extract integration period and time for each baseline
integrations.BL.n = zeros(nBaselines,1);
integrations.BL.integPeriods = zeros(nBaselines,1);
for iBL = 1:nBaselines
    
    integTime   = method.baselines(iBL).IntegTime; % integration time
    integTime   = str2double(integTime);               % integration time, s
    integPeriod = method.baselines(iBL).IntegPeriod; % integ period (string)
    integPeriod = str2double(integPeriod(3:end));      % integ period, ms
    integrations.BL.n(iBL) = integTime/(integPeriod/1e3);
    integrations.BL.integPeriods(iBL) = integPeriod/1e3; % in seconds

    settleTimeiBL = method.baselines(iBL).MagnetSettleTime;    % settle time, ms
    settleTimeiBL = str2double(settleTimeiBL)/1e3;             % settle time, seconds
    settleTime.BL(iBL) = settleTimeiBL; 

end % for iBL

% extract integration period and time for each OP sequence
integrations.OP.n = zeros(nOnPeaks,1);
integrations.OP.integPeriods = zeros(nOnPeaks,1);
for iOP = 1:nOnPeaks

    integTime   = method.onpeaks(iOP).IntegTime; % integration time, s
    integTime   = str2double(integTime);               % integration time, s
    integPeriod = method.onpeaks(iOP).IntegPeriod; % integ period (string)
    integPeriod = str2double(integPeriod(3:end));  % integ period, ms
    integrations.OP.n(iOP) = integTime/(integPeriod/1e3);
    integrations.OP.integPeriods(iOP) = integPeriod/1e3; % in seconds

    settleTimeiOP = method.onpeaks(iOP).MagnetSettleTime;    % settle time, ms
    settleTimeiOP = str2double(settleTimeiOP)/1e3;           % settle time, seconds
    settleTime.OP(iOP) = settleTimeiOP; 

end

nCyclesPerBlock = str2double(method.settings.TotalCycles);
flyBack = method.settings.MagnetFlybackSettleTime; % + time to fly back, ms
settleTime.flyBack = str2double(flyBack)/1e3;

end % function getMethodTiming



