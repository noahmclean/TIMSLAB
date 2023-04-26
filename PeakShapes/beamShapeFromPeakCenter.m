%% script to fit a beam shape and location to a peak center file
% INPUTS 
% - some information about the mass spectrometer e.g. magnet radius
% - location of the .txt file containing peak center output
% 
% OUTPUTS
% - beam shape as a spline
% - beam width/skewness/kurtosis
% - modeled peak shape to compare with measured peak
%
% for Tripoli project by Noah McLean
% 24 April 2023

%% Setup

addpath("exampleData/")
massSpec = massSpecModel("PhoenixPurdue_ATONA");
%filename = "DVCC18-9 z9 Pb-570-PKC-205Pb-PM-S2B7C1.txt";
%filename = "HY30ZK z10 Pb-1004-PKC-205Pb-PM-S2B7C1.txt";
%filename = "HY30ZK z10 Pb-1004-PKC-207Pb-PM-S4B8C1.TXT";
%filename = "6NHCl dpblank-210204A-169-PKC-208Pb-PM-S5B2C1.txt";
%filename = "NBS987 StaticAxH1H2 Bead1Run1-393-PKC-86Sr-Ax-S1B8C1.txt";
%filename = "Ryan-PKC-270UO-PM-Peak.TXT";
filename = "A519_Pb-2531-PKC-208Pb-H2-S3B6C1.TXT";
data = dataModel(filename);


%% This goes sideways for unclear reasons.

%[bestCollectorWidthMM, resnorm] = fminbnd(@(collectorWidthMM) ...
%                                  fitPKCData(massSpec, data, collectorWidthMM), 0, 2);

nWidths = 5000;
resnormVector = zeros(nWidths,1);
widthVector = linspace(0.01, 2, nWidths)';
for iWidth = 1:nWidths
    
    resnormVector(iWidth) = fitPKCData(massSpec, data, widthVector(iWidth));

end
plot(widthVector, resnormVector)


