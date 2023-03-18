%% make synthetic .txt data file for a Phoenix TIMS
% use real method files parsed by FaradayRelativeEfficiency codes
% use a mass spectrometer model from PeakShapes class definition
% add an ion beam intensity time-series
% output: .TXT file in syndata subfolder
% 
% created by Noah McLean for Tripoli and FaradayRelativeEfficiencies 
% on 28-Feb-2023

%% 1. Setup 

% add TIMSLAB to path to use its functions/classes
addpath(genpath("../../../TIMSLAB"));
massSpec = massSpecModel("PhoenixKansas_1e12");

refmat = referenceMaterial("NBS981");
methodName = "Pb 4-5-6-7-8 Daly 10-5-5-5-2 sec.TIMSAM";
nBlocks = 10;
intensityFunction = @(t) 1e6; % cps of major isotope
%intensityFunction = @(ampl, freq, minInt, t)  ...
%      ampl*(ceil(freq*t)-freq*t)+minInt; % sawtooth

% isotopic fractionation for Faradays and Ion Counters
betaFaraday = -0.2; % 0.10%/amu at Pb mass
betaDaly    = -0.3; % 0.15%/amu at Pb mass
% using (a/b)meas = (a/b)true*(Ma/Mb)^beta

% collector relative efficiencies
%          L5  L4  L3  L2  Ax  IC  H1  H2  H3  H4
CREtrue = [1   1   1   1   1   0.9 1   1   1   1];

method = parseTIMSAM(methodName);
CollNames = ["L5", "L4", "L3", "L2", "Ax", "PM", "H1", "H2", "H3", "H4"];
method = processMethod(method, CollNames);

%% 2. Piece out timing from method

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

    integrations.OP

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

end

nCyclesPerBlock = str2double(method.settings(5).Value);
%integrations.OP.settleTime = 
