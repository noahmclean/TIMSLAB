function setup = setupSynDataParams()
%SETUPSYNDATAPARAMS Create setup struct for synthetic data
%   Used to create synthetic .txt file of Phoenix data
%   
%   Created for Tripoli and CollectorRelativeEfficiencies on 5-Apr-2023
%   by Noah McLean


%% Specify synthetic data parameters

% % create a new sample (or use a reference material)
% s.name    = "NBS981Measurement_100kcps206Pb";
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
betaDaly    = @(t) -0.32*ones(size(t)); % 0.16%/amu at Pb mass
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

% write data to file?
writeData = true;

%% Save off timestamp

% log the time the synthetic data is sampled
timeCreated = string(datetime("now", "format", "d MMMM yyyy HH:mm:ss.SSS"));


%% save into struct

setup.sample = s;
setup.spl = spl;
setup.synDataFileName = synDataFileName;
setup.massSpec = massSpec;

setup.nBlocks = nBlocks;
setup.intensityFunction = intensityFunction;
setup.betaFaraday = betaFaraday;
setup.betaDaly = betaDaly;
setup.CREtrue = CREtrue;

setup.darkNoise = darkNoise;
setup.refVolts = refVolts;
setup.kB = kB;
setup.tempInK = tempInK;
setup.method = method;

setup.tStart = tStart;
setup.tBetweenBlocks = tBetweenBlocks;
setup.writeData = writeData;
setup.timeCreated = timeCreated;

end % function setupSynDataParams