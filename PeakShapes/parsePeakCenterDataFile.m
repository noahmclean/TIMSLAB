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
data = readtable(filename, opts);

% Convert to output type
data = table2array(data);
