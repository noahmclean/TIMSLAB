%% make synthetic .txt data file for a Phoenix TIMS
% use real method files parsed by FaradayRelativeEfficiency codes
% use a mass spectrometer model from PeakShapes class definition
% add an ion beam intensity time-series
% output: .TXT file in syndata subfolder
% 
% created by Noah McLean for Tripoli and FaradayRelativeEfficiencies 
% on 28-Feb-2023

%% 1. Setup 

setup = setupSynDataParams();
synDataFileName = setup.synDataFileName;
trueDataFileName = synDataFileName + "_TV";
modelParamFileName = synDataFileName + "_MP";

%% 2. Piece out timing from method

[integrations, settleTime] = getMethodTiming(setup.method);

%% 3. Write the header

writeSynDataHeader(setup, synDataFileName)
writeSynDataHeader(setup, trueDataFileName)

%% 4. Write the collector block

writeSynDataCollector(setup, synDataFileName)
writeSynDataCollector(setup, trueDataFileName)

%% 5. Simulate block data, save off in matrices

idx = makeIndices(setup, integrations, setup.method);
% make a string array for BL and OP data (serial columns and intensities)

BLmeas = [];
BLtrue = [];
OPmeas = [];
OPtrue = [];

tCurrent = setup.tStart;
for iBlock = 1:setup.nBlocks

    % baselines first
    [BLmeas, BLtrue, tCurrent] = makeSyntheticBLData(setup, integrations, ...
                                          settleTime, idx, iBlock, BLmeas, BLtrue, tCurrent);

    % onpeaks next
    [OPmeas, OPtrue, tCurrent] = makeSyntheticOPData(setup, integrations, ...
                                          settleTime, idx, iBlock, OPmeas, OPtrue, tCurrent);

    % allow time for inter-block operations
    tCurrent = tCurrent + setup.tBetweenBlocks;

end % for iBlock = 1:nBlocks


%% 6. Write baselines and onpeaks to file

writeSynDataBLOP(setup, BLmeas, OPmeas, synDataFileName);
writeSynDataBLOP(setup, BLtrue, OPtrue, trueDataFileName);

%% 7. Write true model parameters to file



%% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%

function [integrations, settleTime] = getMethodTiming(method)

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


flyBack = method.settings.MagnetFlybackSettleTime; % + time to fly back, ms
settleTime.flyBack = str2double(flyBack)/1e3;

end % function getMethodTiming


%% MakeIndices

function idx = makeIndices(setup, integrations, method)

nIonCounters = size(setup.massSpec.ionCounterNames,2);
nFaradays = size(setup.massSpec.faradayNames,2);
nCollectors = nIonCounters + nFaradays;
nSerialColumns = 7;
nCyclesPerBlock = str2double(method.settings.TotalCycles);

idx.totalIntegrationsBL = sum(integrations.BL.n)*nCyclesPerBlock*setup.nBlocks;
idx.totalIntegrationsOP = sum(integrations.OP.n)*nCyclesPerBlock*setup.nBlocks;
idx.nCollectors = nIonCounters + nFaradays;
idx.ionCounterColumnIndices = nSerialColumns+1:nSerialColumns+nIonCounters;
idx.faradayColumnIndices = nSerialColumns+nIonCounters+1:nSerialColumns+nCollectors;
idx.ionCounterMethodIndices = 1:nIonCounters; % indices of ion counters in method/F_ind
idx.faradayMethodIndices = nIonCounters + 1: nCollectors; % indices of faradays in F_ind

idx.nIonCounters = nIonCounters;
idx.nFaradays = nFaradays;
idx.nCollectors = nCollectors;
idx.nSerialColumns = nSerialColumns;
idx.nCyclesPerBlock = nCyclesPerBlock;

end % function makeIndices


%% Synthetic Data for Baseline

function [BLmeas, BLtrue, tCurrent] = makeSyntheticBLData(setup, integrations, ...
                                          settleTime, idx, iBlock, BLmeas, BLtrue, tCurrent)

% unpack setup and idx terms
nSerialColumns          = idx.nSerialColumns;
nCollectors             = idx.nCollectors;
ionCounterColumnIndices = idx.ionCounterColumnIndices;
faradayColumnIndices    = idx.faradayColumnIndices;
darkNoise               = setup.darkNoise;

method   = setup.method;
kB       = setup.kB;
tempInK  = setup.tempInK;
massSpec = setup.massSpec;
refVolts = setup.refVolts;

% for each baseline sequence in iBlock
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

    % make separate string arrays for measured and true intensities
    BLseqMeas = BLseq;
    BLseqTrue = BLseq;

    % ion counter baselines (dark noise), units of cps
    lambda = repmat(darkNoise*integrationPeriod, nIntegrations,1);
    darkNoiseIC = random('poisson', lambda);
    BLseqTrue(:,ionCounterColumnIndices) = compose("%1.12e", lambda);
    BLseqMeas(:,ionCounterColumnIndices) = compose("%1.12e", darkNoiseIC);

    % faraday baselines (Johnson noise), units of volts
    s2 = 4 * kB * tempInK * massSpec.amplifierResistance / integrationPeriod;
    mu = refVolts(tvector);
    sigma = repmat(sqrt(s2), nIntegrations, 1);
    johnsonNoiseF = random('normal', mu, sigma);
    BLseqTrue(:,faradayColumnIndices) = compose("%1.12e", mu);
    BLseqMeas(:,faradayColumnIndices) = compose("%1.12e", johnsonNoiseF);

    % update BL with this sequence
    BLmeas = [BLmeas; BLseqMeas]; %#ok<AGROW>
    BLtrue = [BLtrue; BLseqTrue]; %#ok<AGROW>

    % calculate true baseline


end % for iBLseq

end % function makeSyntheticBLData


%% Synthetic Data for OnPeaks
function [OPmeas, OPtrue, tCurrent] = makeSyntheticOPData(setup, integrations, ...
    settleTime, idx, iBlock, OPmeas, OPtrue, tCurrent)

% unpack setup and idx terms
nSerialColumns          = idx.nSerialColumns;
nCollectors             = idx.nCollectors;
ionCounterColumnIndices = idx.ionCounterColumnIndices;
faradayColumnIndices    = idx.faradayColumnIndices;
ionCounterMethodIndices = idx.ionCounterMethodIndices;
faradayMethodIndices    = idx.faradayMethodIndices;

spl               = setup.spl;
massSpec          = setup.massSpec;
method            = setup.method;
kB                = setup.kB;
tempInK           = setup.tempInK;
intensityFunction = setup.intensityFunction;
betaFaraday       = setup.betaFaraday;
betaDaly          = setup.betaDaly;
CREtrue           = setup.CREtrue;
darkNoise         = setup.darkNoise;
refVolts          = setup.refVolts;

% for each OP cycle in iBlock
for iCycle = 1:idx.nCyclesPerBlock
    
    % for each sequence in iCycle
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

        % create OPseq for measured and true values
        OPseqMeas = OPseq;
        OPseqTrue = OPseq;

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
        betaMatrixDaly    = repmat(betaDaly(tvector),    1, idx.nIonCounters);
        betaMatrixFaraday = repmat(betaFaraday(tvector), 1, idx.nFaradays);
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
        %collectEffSeq = -inf(nIntegrations, nCollectors);
        %collectEffSeq(:,collectorRefsInMethod) = CREtrueSeq(:,collectorIndicesUsed);

        % no tails yet, will be a fucntion of 
        tailsSeq = zeros(nIntegrations, nCollectors);

        % 4. put it all together, one step at a time
        % 4a. intensities, isotopically fractionated, + tails
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
            intICEffFractSeq(:,ionCounterMethodIndices) .* dtUnCorrFactor;

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
        OPseqMeas(:,[ionCounterColumnIndices faradayColumnIndices]) = compose("%1.12e", outputOP);

        % d = G(m), calculate 'true' values of measured data
        refCounts = refVolts(tvector) ./ repmat(massSpec.voltsPerCPS, nIntegrations, 1); % cps on Faradays
        OPseqTrueArray = intFractSeq .* CREtrueSeq + [lambda refCounts];  % int * CRE + refs on IC and Faradays
        % convert cps on Faradays to Volts
        OPseqTrueArray(:,faradayMethodIndices) = OPseqTrueArray(:,faradayMethodIndices) .* ...
                                                       repmat(massSpec.voltsPerCPS, nIntegrations, 1);
        % append strings to string array with serial columns
        OPseqTrue(:,[ionCounterColumnIndices faradayColumnIndices]) = compose("%1.12e", OPseqTrueArray);

        % update OP with this sequence
        OPmeas = [OPmeas; OPseqMeas]; %#ok<AGROW>
        OPtrue = [OPtrue; OPseqTrue]; %#ok<AGROW>

    end % for iOPseq

    tCurrent = tCurrent + settleTime.flyBack;

end % for iCycle

end % function makeSyntheticOPData


%% Write synethetic BL and OP to file

function writeSynDataBLOP(setup, BLmeas, OPmeas, filename)

massSpec = setup.massSpec;

%% Write baselines to file

writematrix("#BASELINES", "../syndata/" + filename, ...
    "QuoteStrings", "none", "WriteMode", "append")

writematrix(["ID","Block","Cycle","Integ","PeakID","AxMass","Time",  ...
    massSpec.ionCounterNames, ...
    massSpec.faradayNames], "../syndata/" + filename, ...
    "QuoteStrings", "none", "WriteMode", "append")

% add some spaces after commas in BL
spaceArray = strings(size(BLmeas));
spaceArray(:,2:end) = " ";
BLmeas = spaceArray + BLmeas;

writematrix(BLmeas, "../syndata/" + filename, ...
    "QuoteStrings", "none", "WriteMode", "append")

writematrix(" ", "../syndata/" + filename, ...
    "QuoteStrings", "none", "WriteMode", "append")


%% Write onPeaks to file

writematrix("#ONPEAK", "../syndata/" + filename, ...
    "QuoteStrings", "none", "WriteMode", "append")

writematrix(["ID","Block","Cycle","Integ","PeakID","AxMass","Time",  ...
    massSpec.ionCounterNames, ...
    massSpec.faradayNames], "../syndata/" + filename, ...
    "QuoteStrings", "none", "WriteMode", "append")

% add some spaces after commas in BL
spaceArray = strings(size(OPmeas));
spaceArray(:,2:end) = " ";
OPmeas = spaceArray + OPmeas;

writematrix(OPmeas, "../syndata/" + filename, ...
    "QuoteStrings", "none", "WriteMode", "append")

writematrix([" "; "#END,AnalysisCompleted"], "../syndata/" + filename, ...
    "QuoteStrings", "none", "WriteMode", "append")

end % function writeSynDataOPBL