classdef dataModel
    %DATAMODEL Model of peak center/centre data file information
    %   Data includes measured masses and intensities
    %   created by Noah McLean 25 April 2022
    %   for Peak Center/Beam Shape project
    
    properties
        magnetMasses      % vector of masses for intensity measurements
        measPeakIntensity % vector of corresponding peak intensities
        peakCenterMass    % mass at center of peak from header
        integPeriodMS     % integration period of measurements in ms
        MassID            % name of peak getting centered e.g. "205Pb"
        detectorName      % name of detector as string e.g. "L2"
        collectorWidthAMU % width of collector aperture in AMU at center mass
        beamWidthAMU      % width of beam in AMU at center mass

    end
    
    methods
        function data = dataModel(filename)
            %DATAMODEL Construct an instance of this class
            %   Detailed explanation goes here

            arguments 
                filename (1,1) string
            end 

            dataTable = parseDataFromTextFile(filename);
            data.magnetMasses = dataTable.Mass;
            data.measPeakIntensity = dataTable.Intensity;

            % 2. parse header and filename
            
            % read header as string array
            fileAsStrings = readlines(filename, 'EmptyLineRule', 'skip');
            header = fileAsStrings(1:11); % just the header
            
            % extract useful parts of header
            data.peakCenterMass = str2double(extractAfter(header(5), ","));
            data.integPeriodMS = str2double(extractAfter(header(11), "ms"));
            data.MassID = strtrim(extractAfter(header(3), ","));
            
            % read detector name from filename, could do from header too?
            pat = regexpPattern('\w*'); % regular expression for 'words'
            filenameBits = extract(filename, pat); % etract words from pattern
            data.detectorName = filenameBits(end-2); % second to last word

        end % cunstructor function

        function collectorWidthAMU = calcCollectorWidthAMU(data, massSpec)
            %CALCCOLLECTORWIDTHAMU Collector width in AMU
            %   Calculate collector aperture width in AMU at measured mass
            %   massAtCenterAMU is average/peak mass of the scan, in amu
            
            collectorWidthAMU = data.peakCenterMass / ...
                massSpec.effectiveRadiusMagnetMM * massSpec.collectorWidthMM;

        end % calcCollectorWidthAMU

        function beamWidthAMU = calcBeamWidthAMU(data, massSpec)
            %CALCBEAMWIDTHAMU Beam width in AMU
            %   Calculate beam width in AMU at measured mass
            %   massAtCenterAMU is average/peak mass of the scan, in amu

            beamWidthAMU = data.peakCenterMass / ...
                massSpec.effectiveRadiusMagnetMM * massSpec.theoreticalBeamWidthMM;
        end % calcBeamWidthAMU
        
    end % methods

end % classdef


function dataTable = parseDataFromTextFile(filename)
% PARSEDATAFROMTEXTFILE Parse text file and return only data table

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

end % parse text file and
