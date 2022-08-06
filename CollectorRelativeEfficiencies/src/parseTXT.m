function data = parseTXT(dataFolder)
%PARSETXT Parse Isotopx TXT file output
%   dataFolder is .RAW folder, method struct from parseTIMSAM

textFileInfo = dir("../data/" + dataFolder + "/*.TXT");

%% parse header block

% get method row from header, find method name
opts = delimitedTextImportOptions('NumVariables', 2);
opts.DataLines = [4, 11];
methodHeader = readcell(textFileInfo.name, opts);

data.header.fileName   = string(methodHeader{1,2});
data.header.methodName = string(methodHeader{2,2});
data.header.methodPath = string(methodHeader{3,2});
data.header.IsoWorksMethod = string(methodHeader{4,2});
data.header.FolderPath = string(methodHeader{5,2});
data.header.Corrected = string(methodHeader{6,2});
data.header.BChannels = string(methodHeader{7,2});
data.header.TimeZero = string(methodHeader{8,2});


%% new "parse data" required for ATONA A/B channel data

opts = delimitedTextImportOptions;
opts.VariableTypes = "string";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
firstColumn = readmatrix(textFileInfo.name, opts);
% note that default is to skip empty lines
collectorsStartPosition = find(firstColumn == "#COLLECTORS");
userTablesStartPosition = find(firstColumn == "#USERTABLES");
baselinesStartPosition  = find(firstColumn == "#BASELINES");
onPeakStartPosition     = find(firstColumn == "#ONPEAK");
endPosition             = find(firstColumn == "#END");

collectorsEndPosition = min([userTablesStartPosition baselinesStartPosition onPeakStartPosition]) - 2;
nCollectors = collectorsEndPosition - collectorsStartPosition - 1;
data.collectorNames = firstColumn(collectorsStartPosition+2:collectorsEndPosition);

if data.header.BChannels == "No" % if resistor-based amplifiers, no BChannels
    nDataColumns = 7 + nCollectors;
elseif data.header.BChannels == "Yes" % if ATONAs
    nDataColumns = 7 + 2*nCollectors - 1; % minus one because PM doesn't have BChannel
else
    disp('unrecognized text file column setup')
end

% extract collector block information
opts.VariableNames = string(1:6);
opts.VariableTypes = ["string", "string", "string", "string", "string", "string"];
opts.DataLines = [collectorsStartPosition+1 collectorsEndPosition]; 
data.Collectors = readmatrix(textFileInfo.name, opts);
% including gains and resistances
firstFaradayRow = find(data.Collectors(:,2) == "F", 1, "first");
 lastFaradayRow = find(data.Collectors(:,2) == "F", 1, "last");
 Frange = firstFaradayRow:lastFaradayRow;
data.FaradayResist = double(data.Collectors(Frange,3)); % resistances (ohms)
data.FaradayGains = double(data.Collectors(Frange,4));  % gains (relative to Axial)

% next: line 104 (parse BL and OP data, separate out indices and data)

% baselines
if ~isempty(baselinesStartPosition) % if baselines exist
opts.VariableNames = string(1:nDataColumns);
opts.VariableTypes = repmat("string", 1, nDataColumns);
opts.DataLines = [baselinesStartPosition+2 onPeakStartPosition-2]; 
data.BLall = readmatrix(textFileInfo.name, opts);

data.BLserial = double(data.BLall(:,2:4));% [block cycle integration] serially assigned counts

data.BLmatrix = double(data.BLall(:,8:end)); % matrix of collector readings
if data.header.BChannels == "Yes"
    data.BLmatrixA = data.BLmatrix(:,1:nCollectors);
    data.BLmatrixB = data.BLmatrix(:,nCollectors+1:end);
end

data.BLtime   = double(data.BLall(:,7)); % time
data.BLID = data.BLall(:,1); % baseline ID, eg "BL1", "BL2", etc. 1st column in TXT data file
data.BLSeqIdx = double(extractAfter(data.BLall(:,1), "BL"));
end % if baselines exist

% on-peaks
opts.DataLines = [onPeakStartPosition+2 endPosition-2]; 
data.OPall = readmatrix(textFileInfo.name, opts);

data.OPserial = double(data.OPall(:,2:4));% [block cycle integration] serially assigned counts

data.OPmatrix = double(data.OPall(:,8:end)); % matrix of all collector readings
if data.header.BChannels == "Yes" % if ATONAs
    data.OPmatrixA = data.OPmatrix(:,1:nCollectors);
    data.OPmatrixB = data.OPmatrix(:,nCollectors+1:end);
end

data.OPtime   = double(data.OPall(:,7)); % time
data.OPID = data.OPall(:,1); % OnPeak ID, eg "OP1", "OP2", etc.  1st column in TXT data file
data.OPSeqIdx = double(extractAfter(data.OPall(:,1), "OP"));

% 
% %% parse data 
% 
% % future work: change to table import for ATONA data
% %opts = setvartype(opts, 'string');
% %opts.DataLines = [13 Inf];
% %dataTable = readtable(textFileInfo.name, opts);
% 
% fid=fopen(textFileInfo.name,'r');
% dtmp=textscan(fid,'%s','delimiter',',','Headerlines',12); 
% dtmp = string(dtmp{1});
% fclose(fid);
% 
% collectorsStartPosition = find(dtmp == "#COLLECTORS");
% userTablesStartPosition = find(dtmp == "#USERTABLES");
% baselinesStartPosition  = find(dtmp == "#BASELINES");
% onPeakStartPosition     = find(dtmp == "#ONPEAK");
% endPosition             = find(dtmp == "#END");
% 
% % determine number of data columns by counting collectors
% % there are 6 columns in COLLECTORS table, 1 row of header
% collectorBlockEndPosition = min([userTablesStartPosition baselinesStartPosition onPeakStartPosition]) - 1;
% nCollectors = (collectorBlockEndPosition-collectorsStartPosition)/6 - 1;
% data.collectorNames = dtmp(8:6:collectorBlockEndPosition)'; % starts at 
% 
% if data.header.BChannels == "No" % if resistor-based amplifiers, no BChannels
%     nDataColumns = 7 + nCollectors;
% elseif data.header.BChannels == "Yes" % if ATONAs
%     nDataColumns = 7 + 2*nCollectors - 1; % minus one because PM doesn't have BChannel
% else
%     disp('unrecognized text file column setup')
% end
% 
% % grab the gains
% collRange = (collectorsStartPosition+1):collectorBlockEndPosition;
% data.Collectors = reshape(dtmp(collRange), 6, [])';
% firstFaradayRow = find(data.Collectors(:,2) == "F", 1, "first");
%  lastFaradayRow = find(data.Collectors(:,2) == "F", 1, "last");
%  Frange = firstFaradayRow:lastFaradayRow;
% data.FaradayResist = double(data.Collectors(Frange,3)); % resistances (ohms)
% data.FaradayGains = double(data.Collectors(Frange,4));  % gains (relative to Axial)
% 
% % range starts after header, continues to cell before next block flag
% BLrange = (baselinesStartPosition+1+nDataColumns):(onPeakStartPosition-1);
% OPrange = (onPeakStartPosition+1+nDataColumns):(endPosition-1);
% 
% data.BLall = reshape(dtmp(BLrange),nDataColumns,[])';
% data.OPall = reshape(dtmp(OPrange),nDataColumns,[])';
% 
% data.BLserial = double(data.BLall(:,2:4));% [block cycle integration] serially assigned counts
% data.BLmatrix = double(data.BLall(:,8:end)); % matrix of collector readings
% data.BLtime   = double(data.BLall(:,7)); % time
% data.BLID = data.BLall(:,1); % baseline ID, eg "BL1", "BL2", etc. 1st column in TXT data file
% data.BLSeqIdx = double(extractAfter(data.BLall(:,1), "BL"));
% 
% data.OPserial = double(data.OPall(:,2:4));% [block cycle integration] serially assigned counts
% data.OPmatrix = double(data.OPall(:,8:end)); % matrix of collector readings
% data.OPtime   = double(data.OPall(:,7)); % time
% data.OPID = data.OPall(:,1); % OnPeak ID, eg "OP1", "OP2", etc.  1st column in TXT data file
% data.OPSeqIdx = double(extractAfter(data.OPall(:,1), "OP"));
% 




end % function parseTXT

