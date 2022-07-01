function tails = initializePeakTails(method)
%INITIALIZEPEAKTAILS Summary of this function goes here
%   Taken from peakTailsVisualization_v2.m by Noah McLean for Tripoli
%   measurement model specification.
%   tails contains peak tail intensities as a proportion of the most
%   abundant isotope measured in the input argument method
%   tails.OP is expected on-peak tails contribution as a proportion of
%   the major isotope from the method.  
%   tails.halfMass is half-mass tails contributions at 0.5 amu above
%   and below the tails.OP peaks.

element = extract( string(method.onpeaks(1).Info(3).Value), lettersPattern);
switch element
    case "Sm"
    isotopes = ["144Sm", "147Sm", "148Sm", "149Sm", "150Sm", "152Sm", "154Sm"];
    abundances = [3.08 15 11.25 13.82 7.37 26.74 22.74];
    case "Pb"
    isotopes = ["204Pb", "206Pb", "207Pb", "208Pb"];
    abundanceRatios = [0.0590075 1 0.914583 2.1681]; % Condon
    abundances = abundanceRatios/sum(abundanceRatios)*100;
end % switch


% key to mass.m value class with isotopic masses
massNames = element + extractBefore(isotopes, element);
peaks = zeros(size(isotopes));
for iPeak = 1:length(isotopes)
    peaks(iPeak) = mass.(massNames(iPeak));
end

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

% get tails for the method's isotopes. extract method info
massIDs = method.MassIDs';
nMassIDs = length(massIDs);
massIDNames = element + extractBefore(massIDs, element);

% determine relative abundances of method's massIDs to renormalize if needed
methodIntensity = zeros(1,nMassIDs);
for iMass = 1:nMassIDs
    methodIntensity(iMass) = intensities(massNames == massIDNames(iMass));
end % for iMass determine relative abundances
[maxMethodIntensity, maxMethodIntensityIdx] = max(methodIntensity);
%rescale sumTails to maxMethodIntensity, from maxIntensity for all isotopes
sumTails = sumTails * maxIntensity/maxMethodIntensity;

% interpolate sumTails at isotopic masses (OP) and half-masses
% tails is 
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

