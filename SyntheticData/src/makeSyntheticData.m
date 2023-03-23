%% make synthetic .txt data file for a Phoenix TIMS
% use real method files parsed by FaradayRelativeEfficiency codes
% use a mass spectrometer model from PeakShapes class definition
% add an ion beam intensity time-series
% output: .TXT file in syndata subfolder
% 
% created by Noah McLean for Tripoli and FaradayRelativeEfficiencies 
% on 28-Feb-2023

%% 1. Setup 

s.name    = "NBS981Measurement";
s.element = "Pb";
s.species =            ["204Pb",  "205Pb", "206Pb", "207Pb",  "208Pb"];
s.relativeAbundances = [0.0590074, 1e-6,   1,       0.914683, 2.1681];
spl = sample(s.name, s.element, s.species, s.relativeAbundances);

% name the data file -- refactor?
synDataFileName = s.name;

% add TIMSLAB to path to use its functions/classes
addpath(genpath("../../../TIMSLAB"));
massSpec = massSpecModel("PhoenixKansas_1e12");
methodName = "Pb 4-5-6-7-8 Daly 10-5-5-5-2 sec.TIMSAM";
nBlocks = 10;

intensityFunction = @(t) 1e6*ones(size(t)); % cps of major isotope
%intensityFunction = @(ampl, freq, minInt, t)  ...
%      ampl*(ceil(freq*t)-freq*t)+minInt; % sawtooth

% isotopic fractionation for Faradays and Ion Counters
betaFaraday = @(t) -0.2; % 0.10%/amu at Pb mass
betaDaly    = @(t) -0.3; % 0.15%/amu at Pb mass
% using (a/b)meas = (a/b)true*(Ma/Mb)^beta

% collector relative efficiencies
%          PM  RS L5  L4  L3  L2  Ax  H1  H2  H3  H4
CREtrue = [0.9 1  1   1   1   1   1   1   1   1   1];

% baselines
darkNoise = [0.2 0]; % cps, [PM RS]
%refVolts:   L5    L4    L3    L2    Ax     H1    H2    H3    H4
refVolts  = [-1e-2 -2e-2 -1e-2 -2e-2 -3e-2 -1e-2 -2e-2 -1e-2 -2e-2];
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
%BL = strings(totalIntegrationsBL, nSerialColumns+nCollectors);
%OP = strings(totalIntegrationsOP, nSerialColumns+nCollectors);
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
        noisedata = random('poisson', lambda);
        BLseq(:,ionCounterColumnIndices) = compose("%1.12e", noisedata);

        % faraday baselines (Johnson noise), units of volts
        s2 = 4 * kB * tempInK * massSpec.amplifierResistance / integrationPeriod;
        mu = repmat(refVolts, nIntegrations, 1);
        sigma = repmat(sqrt(s2), nIntegrations, 1);
        noisedata = random('normal', mu, sigma);
        BLseq(:,faradayColumnIndices) = compose("%1.12e", noisedata);

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
        OPseq(:,5) = string(method.onpeaks(iOPseq).Info(3).Value); % PeakID
        OPseq(:,6) =  num2str(method.axialMasses.OP(iOPseq), '%1.4f');

        % times for each integration; update tCurrent
        tCurrent = tCurrent + settleTime.OP(iOPseq); % settle before 1st integration
        tStop = tCurrent + (nIntegrations-1)*integrationPeriod;
        tvector = linspace(tCurrent, tStop, nIntegrations)';
        OPseq(:,7) = num2str(tvector, '%1.7f');
        tCurrent = tStop + integrationPeriod;

        % on peaks: ion counters
        % 1. dark noise
        lambda = repmat(darkNoise*integrationPeriod, nIntegrations,1);
        noisedata = random('poisson', lambda);
        % 2. scaled intensities
        speciesIntensities = intensityFunction(tvector) * spl.relativeAbundances;
        % 3. sequence map
        collectorRefsForSequence = method.F_ind(iOPseq,ionCounterMethodIndices);
        % RESUME HERE
        %collectorIdcsForSequence = intensityFunction(tvector) .* 
        
        % x = [1 3 5 7 9]                 vector of scaled intensities
        % y = [0 0 0 0 2 0 0 3 0 0 4 0 5] F_ind
        % z = zeros(1,13)
        % z(find(y>0)) = x(y(y>0)) % distribute elements of x to positions in y

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
    
    integTime   = method.baselines(iBL).Info(6).Value; % integration time
    integTime   = str2double(integTime);               % integration time, s
    integPeriod = method.baselines(iBL).Info(5).Value; % integ period (string)
    integPeriod = str2double(integPeriod(3:end));      % integ period, ms
    integrations.BL.n(iBL) = integTime/(integPeriod/1e3);
    integrations.BL.integPeriods(iBL) = integPeriod/1e3; % in seconds

    settleTimeiBL = method.baselines(iBL).Info(9).Value;    % settle time
    settleTimeiBL = str2double(settleTimeiBL)/1e3;             % settle time, seconds
    settleTime.BL(iBL) = settleTimeiBL; 

end % for iBL

% extract integration period and time for each OP sequence
integrations.OP.n = zeros(nOnPeaks,1);
integrations.OP.integPeriods = zeros(nOnPeaks,1);
for iOP = 1:nOnPeaks

    integTime   = method.onpeaks(iOP).Info(7).Value; % integration time, s
    integTime   = str2double(integTime);               % integration time, s
    integPeriod = method.onpeaks(iOP).Info(6).Value; % integ period (string)
    integPeriod = str2double(integPeriod(3:end));  % integ period, ms
    integrations.OP.n(iOP) = integTime/(integPeriod/1e3);
    integrations.OP.integPeriods(iOP) = integPeriod/1e3; % in seconds

    settleTimeiOP = method.onpeaks(iOP).Info(14).Value;    % settle time
    settleTimeiOP = str2double(settleTimeiOP)/1e3;             % settle time, seconds
    settleTime.OP(iOP) = settleTimeiOP; 

end

nCyclesPerBlock = str2double(method.settings(5).Value);
flyBack = method.settings(13).Value;
settleTime.flyBack = str2double(flyBack)/1e3;

end % function getMethodTiming



