%% Illustration of peak tail corrections.

%% Input 

% define an element of interest
% element = "Pb";
% peaks = [204 206 207 208]; % integer isotopic masses
% abundances = [1.4 24.1 22.1 52.4]/100; % as fraction
% element = "Nd";
% peaks = [142 143 144 145 146 148 150];
% abundances = [27.2 12.2 23.8 8.3 17.2 5.8 5.6]/100;
element = "Sm";
peaks = [144 147 148 149 150 152 154];
abundances = [3.08 15 11.25 13.82 7.37 26.74 22.74];
% element = "U";
% peaks = [233 234 235 236 238];
% abundances = [1 0.005 0.72 1 99.274]/100;

maxTailIntensity = 1.5e-5;
halfPeakWidth = 0.3;

% peak tail parameters
downMassTail = @(mass, peakMass, peakIntensity, a, b)  a * (peakMass-mass).^b * peakIntensity;
upMassTail   = @(mass, peakMass, peakIntensity, c, d)  c * (mass-peakMass).^d * peakIntensity;

% set up parameters of analysis/system
nMasses = 1000;
massRange = linspace(min(peaks)-1.5, max(peaks)+1.5, nMasses);
maxIntensity = 1; % Volts or other unit, for highest abundance isotope

% calculate intensities of peaks
intensities = abundances * maxIntensity/max(abundances);

% reasonable constants for Faradays that are not behind a filter
a = 2e-6; % 2 ppm down-mass abundance sensitivity
b = -2;   % slope governs how fast the peak tail tapers off
c = 0.6e-6; % 1 ppm up-mass abundance sensitivity
d = -3;   % up-mass peak tail tapers more quickly.

printWithSum = true;

%% Calculate peak tail intensities

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

%% Plot peak tails

% figure and axes creation
figureHandle = figure('Position', [1 1 1000 700]);
axLoInt = axes(figureHandle, 'Units', 'normalized', 'Position', [0.1 0.11 0.8 0.55]);
axHiInt = axes(figureHandle, 'Units', 'normalized', 'Position', [0.1 0.70 0.8 0.22]);
hold([axLoInt, axHiInt], 'on')

% axes setup - labels, font sizes, etc.
tickLabelFontSize = 20;
axisLabelFontSize = 24;
lineWidth = 3;
set(axHiInt, 'XTickLabels', {}, 'FontSize', tickLabelFontSize)
set(axLoInt, 'FontSize', tickLabelFontSize, 'YLim', [0 maxTailIntensity])
set([axHiInt, axLoInt], 'XLim', [min(massRange), max(massRange)])
axLoInt.YAxis.Exponent = -6;
xlabel(axLoInt, 'Mass', 'FontSize', axisLabelFontSize)
text(-0.08, 0.7, 'Relative Intensity', ...
                 'FontSize', axisLabelFontSize, ...
                 'Units', 'normalized', 'Rotation', 90, 'Parent', axLoInt)
axHiInt.Title.String = "Peak Tail Model for " + element;
axHiInt.Title.FontWeight = 'Normal';
axHiInt.Title.FontSize = 28;
axLoInt.LineWidth = 2;
axHiInt.LineWidth = 2;


cmap = lines(nPeaks);

% plot up- and down-mass tails in low-intensity (lower) axes
for iPeak = 1:nPeaks
    
    % down-mass tail for iPeak
    tailRange = massRange < ( peaks(iPeak) - halfPeakWidth );
    plot(axLoInt, massRange(tailRange), downTails(tailRange, iPeak), ...
         'LineWidth', lineWidth, 'Color', cmap(iPeak,:))

    % up-mass tail for iPeak
    tailRange = massRange > ( peaks(iPeak) + halfPeakWidth );
    plot(axLoInt, massRange(tailRange), upTails(tailRange, iPeak), ...
         'LineWidth', lineWidth, 'Color', cmap(iPeak,:))

end % for iPeak

% plot sum, remove vertical breaks
breaks = sort([peaks - halfPeakWidth, peaks + halfPeakWidth]);
breaks = [min(massRange) breaks max(massRange)]; % add min/max mass to segments
for iSegment = 1:length(breaks)-1
    segmentStart = breaks(iSegment);
    segmentEnd   = breaks(iSegment + 1);
    inSegment = massRange >= segmentStart & massRange <= segmentEnd;
    
    if printWithSum
    plot(axLoInt, massRange(inSegment), sumTails(inSegment), ...
              'LineWidth', lineWidth+1, 'Color', 'k')
    end

end % for iSegment


%% Make mass spec peaks for higher-intensity axes
% make the shape of a mass spec peak, centered at x = 0

slitWidth = 0.6;
peakStd = 0.03;
beamMass = linspace(-0.2,0.2,200);
beamInt = normpdf(beamMass, 0, peakStd);
slitInt = double( abs(beamMass)<slitWidth/2);
peakShape = conv(beamInt, slitInt, 'full');
peakShape = peakShape / max(peakShape); % normalize
peakMass = linspace(-0.4, 0.4, 399);
%plot(peakMass, peakShape);

% plot a peak shape for each peak, scaled by its intensity
for iPeak = 1:nPeaks

    peakStart = peaks(iPeak) - halfPeakWidth;
    peakEnd   = peaks(iPeak) + halfPeakWidth;
    iPeakIsInMassRange = massRange >= peakStart & massRange <= peakEnd;
    peakScaled = intensities(iPeak) * interp1(peakMass+peaks(iPeak), peakShape, massRange(iPeakIsInMassRange));
    plot(axHiInt, massRange(iPeakIsInMassRange), peakScaled, 'LineWidth', lineWidth, 'Color', cmap(iPeak,:))

end % for iPeak, plot a shape for each peak

% plot a point at half and integer masses, if printing with sum
if printWithSum
    for iPeak = 1:nPeaks
        plot(axLoInt, peaks(iPeak), interp1(massRange, sumTails, peaks(iPeak)), 'o', ...
            'MarkerFaceColor', cmap(iPeak,:), 'MarkerSize', 12, ...
            'MarkerEdgeColor', 'k')
    end

    halfMassStart = bitor(floor(2*massRange(1)   + 1), 1)/2; % nearest 0.5 above min
    halfMassEnd   = bitor(floor(2*massRange(end) - 1), 1)/2; % nearest 0.5 below max
    halfMasses = halfMassStart:halfMassEnd;
    plot(axLoInt, halfMasses, interp1(massRange, sumTails, halfMasses), 'o', ...
        'MarkerFaceColor', 'k', 'MarkerSize', 12, ...
        'MarkerEdgeColor', 'k')
end % if printWithSum

%% Save as pdf

set(figureHandle,'Units','inches', 'InvertHardCopy', 'off');
screenposition = get(figureHandle,'Position');
set(figureHandle,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
if printWithSum
    print -dpdf -vector peakTailsWithSum
    print -dpng peakTailsWithSum
else
    print -dpdf -vector peakTailsWithoutSum
    print -dpng peakTailsWithoutSum
end
