function data = parsePeakCenterDataFile(filename)
% parse a mass spectrometer file
%
% filename is a string containing the path/name of peak center data file
% data is a two-column matrix contining mass and intensity
%
% created by Noah McLean, 6-April-2022

% set options, parse data file
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [14, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Mass", "Intensity"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
dataTable = readtable(filename, opts);

% Convert to output type
data.meas = table2array(dataTable);


%% parse header

fileAsStrings = readlines(filename, 'EmptyLineRule', 'skip');
header = fileAsStrings(1:11); % just the header

data.peakCenterMass = str2double(extractAfter(header(5), ","));
data.integPeriodMS = str2double(extractAfter(header(11), "ms"));
data.MassID = strtrim(extractAfter(header(3), ","));

filenameBits = extractBetween(filename, "-", "-");
data.detectorName = filenameBits(end); clear filenameBits


end % function