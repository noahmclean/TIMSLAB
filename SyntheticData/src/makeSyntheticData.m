%% make synthetic .txt data file for a Phoenix TIMS
% use real method files parsed by FaradayRelativeEfficiency codes
% use a mass spectrometer model from PeakShapes class definition
% add an ion beam intensity time-series
% output: .TXT file in syndata subfolder
% 
% created by Noah McLean for Tripoli and FaradayRelativeEfficiencies on 28-Feb-2023

%% 1. Setup 

% add TIMSLAB to path to use its functions/classes
addpath(genpath("../../../TIMSLAB"));
massSpec = massSpecModel("PhoenixKansas_1e12");

refmat = referenceMaterial("NBS981");
methodName = "Pb 4-5-6-7-8 Daly 10-5-5-5-2 sec.TIMSAM";
nBlocks = 10;
intensityFunction = @(t) 1e6; % cps of major isotope
%intensityFunction = @(ampl, freq, minInt, t)  ampl*(ceil(freq*t)-freq*t)+minInt; % sawtooth

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

%% 2.  