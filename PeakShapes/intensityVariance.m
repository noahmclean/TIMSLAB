function varInt = intensityVariance(measIntensity, collectorName, massSpec, integrationTimeMS)
%INTENSITYVARIANCE Estimate variance or covariance matrix of measured intensities
%   Use shot noise + (Johnson noise or ATONA noise model)
%   measIntensity is scalar or vector. If vector, return cov. matrix
%   collectorName is the name of the collector (string)
%   massSpec is the massSpecModel being used
%   integrationTimeMS is the integration time in milliseconds
%
%   Assumption: if the collector is a Faraday, its output will be in V
%   (10^11 - equivalent)
%
%   by Noah McLean for the Tripoli and FaradayRelativeEfficiency projects
%   25 April 2023

arguments % validation
    measIntensity      (1,:) double
    collectorName      (1,1) string
    massSpec           (1,1) massSpecModel
    integrationTimeMS  (1,1) double
end

%% 1. identify Faraday vs. Ion Counter

faradayNameMatches    = matches(massSpec.faradayNames,    collectorName);
ionCounterNameMatches = matches(massSpec.ionCounterNames, collectorName);

if any(faradayNameMatches)
    collectorType = "Faraday";
    collectorIndex = find(faradayNameMatches);
elseif any(ionCounterNameMatches)
    collectorType = "ionCounter";
    collectorIndex = find(ionCounterNameMatches);
else
    error("Collector name not recognized")
end % if block to identify collectorType

% initialize variance/covariance matrix for intensity
nData = length(measIntensity);
varInt = zeros(nData);

integrationTimeSeconds = integrationTimeMS/1000;

%% 2. shot noise

if collectorType == "Faraday" % then convert to cps
    cpsFactorsAll = cpsPerVolt(massSpec);
    cpsFactorThisCup = cpsFactorsAll(collectorIndex);
    voltFactorsAll = voltsPerCPS(massSpec);
    voltFactorThisCup = voltFactorsAll(collectorIndex);
else
    cpsFactorThisCup = 1;
end % if a conversion factor needed for Faraday

% measured ion arrivals (ignoring dead time correction)
measIntCounts = measIntensity * cpsFactorThisCup * integrationTimeSeconds;

% shote noise (using 1 as a minimum to keep matrix nonsingular)
shotNoiseCounts = diag(max(1, measIntCounts)); % variance = counts (Poisson distributed)
shotNoiseCPS = shotNoiseCounts * (1/integrationTimeSeconds)^2; % cps = counts/seconds, error propagation

if collectorType == "Faraday" % then convert cps back to volts
    shotNoise = shotNoiseCPS * voltFactorThisCup^2; % squared because error propagation
else 
    shotNoise = shotNoiseCPS; % if ion counter, data are all in cps
end 

varInt = varInt + shotNoise;


%% 3. amplifier noise (if Faraday)

if collectorType == "Faraday"

    switch massSpec.faradayTypes(collectorIndex)
        
        case "resistance"
        R = massSpec.amplifierResistance(collectorIndex);
        ampNoise = 4 * massSpec.kB * massSpec.tempInK * R / integrationTimeSeconds;

        case "ATONA"
        ampNoiseCPS = massSpec.noiseConstantATONAsCPS * 1/integrationTimeSeconds;
        ampNoise = ampNoiseCPS * voltFactorThisCup; % convert cps to volts

    end % switch case restistance vs. ATONA

ampNoise = diag(ampNoise*ones(nData,1));
varInt = varInt + ampNoise;

end % function intensityVariance

