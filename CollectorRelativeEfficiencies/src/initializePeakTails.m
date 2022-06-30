function tails = initializePeakTails(method)
%INITIALIZEPEAKTAILS Summary of this function goes here
%   Taken from peakTailsVisualization_v2.m by Noah McLean for Tripoli
%   measurement model specification

element = "Sm";
isotopes = ["144Sm", "147Sm", "148Sm", "149Sm", "150Sm", "152Sm", "154Sm"];

% key to mass.m value class with isotopic masses
massNames = element + extractBefore(isotopes, element);
peaks = zeros(size(isotopes));
for iPeak = 1:length(isotopes)
    peaks(iPeak) = mass.(massNames(iPeak));
end

abundances = [3.08 15 11.25 13.82 7.37 26.74 22.74];
halfPeakWidth = 0.3;

downMassTail = @(mass, peakMass, peakIntensity, a, b)  a * (peakMass-mass).^b * peakIntensity;
upMassTail   = @(mass, peakMass, peakIntensity, c, d)  c * (mass-peakMass).^d * peakIntensity;

% set up parameters of analysis/system
nMasses = 1000;
massRange = linspace(min(peaks)-1.5, max(peaks)+1.5, nMasses);
maxIntensity = 1; % Volts or other unit, for highest abundance isotope

% normalize intensities of peaks
intensities = abundances * maxIntensity/max(abundances);

% reasonable constants for Faradays that are not behind a filter
a = 2e-6; % 2 ppm down-mass abundance sensitivity
b = -2;   % slope governs how fast the peak tail tapers off
c = 0.6e-6; % 1 ppm up-mass abundance sensitivity
d = -3;   % up-mass peak tail tapers more quickly.

nPeaks = length(peaks);
downTails = zeros(length(massRange), nPeaks);
upTails   = zeros(length(massRange), nPeaks);

for iPeak = 1:nPeaks
    
    % down-mass tail

    % mass range in which tail is defined:
    tailRange                  = massRange < ( peaks(iPeak) - halfPeakWidth );
    % intensity of tail across all masses (will pick out defined range next)
    downTailiMass              = downMassTail(massRange, peaks(iPeak), intensities(iPeak), a, b);
    % save off the tail intensities in the correct rainge = tailRange
    downTails(tailRange,iPeak) = downTailiMass(tailRange);

    % up-mass tail

    % mass range in which tail is defined:
    tailRange                  = massRange > ( peaks(iPeak) + halfPeakWidth );
    % intensity of tail across all masses (will pick out defined range next)
    upTailiMass                = upMassTail(massRange, peaks(iPeak), intensities(iPeak), c, d);
    % save off the tail intensities in the correct rainge = tailRange
    upTails(tailRange,iPeak) = upTailiMass(tailRange);

end % for iPeak

sumTails = sum([upTails downTails], 2);

% get tails for the method's isotopes
massIDs = method.MassIDs';
nMassIDs = length(massIDs);
massIDNames = element + extractBefore(isotopes, element);

tails.OP = zeros(2,nMassIDs);
tails.halfMass = zeros(2,nMassIDs+1);
for iMass = 1:nMassIDs
    
    thisMass = mass.(massIDNames(iMass));
    tails.OP(1, iMass) = thisMass;
    tails.OP(2, iMass) = interp1(massRange, sumTails, thisMass);

end % for iMass

halfMasses = [tails.OP(1,:)-0.5 tails.OP(1,end)+0.5]; % all the half-masses currently used
tails.halfMass = [halfMasses; interp1(massRange, sumTails, halfMasses)];


end

