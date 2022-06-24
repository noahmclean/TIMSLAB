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


%% parse data 

fid=fopen(textFileInfo.name,'r');
dtmp=textscan(fid,'%s','delimiter',',','Headerlines',12); 
dtmp = string(dtmp{1});
fclose(fid);

collectorsStartPosition = find(dtmp == "#COLLECTORS");
userTablesStartPosition = find(dtmp == "#USERTABLES");
baselinesStartPosition  = find(dtmp == "#BASELINES");
onPeakStartPosition     = find(dtmp == "#ONPEAK");
endPosition             = find(dtmp == "#END");

% determine number of data columns by counting collectors
% there are 6 columns in COLLECTORS table, 1 row of header
collectorBlockEndPosition = min([userTablesStartPosition baselinesStartPosition onPeakStartPosition]) - 1;
nCollectors = (collectorBlockEndPosition-collectorsStartPosition)/6 - 1;
data.collectorNames = dtmp(8:6:collectorBlockEndPosition)'; % starts at 

if data.header.BChannels == "No" % if resistor-based amplifiers, no BChannels
    nDataColumns = 7 + nCollectors;
elseif data.header.BChannels == "Yes" % if ATONAs
    nDataColumns = 7 + 2*nCollectors;
else
    disp('unrecognized text file column setup')
end

% range starts after header, continues to cell before next block flag
BLrange = (baselinesStartPosition+1+nDataColumns):(onPeakStartPosition-1);
OPrange = (onPeakStartPosition+1+nDataColumns):(endPosition-1);

data.BLall = reshape(dtmp(BLrange),nDataColumns,[])';
data.OPall = reshape(dtmp(OPrange),nDataColumns,[])';

data.BLserial = double(data.BLall(:,2:4));% [block cycle integration] serially assigned counts
data.BLmatrix = double(data.BLall(:,8:end)); % matrix of collector readings
data.BLtime   = double(data.BLall(:,7)); % time
data.BLID = data.BLall(:,1); % baseline ID, eg "BL1", "BL2", etc.

data.OPserial = double(data.OPall(:,2:4));% [block cycle integration] serially assigned counts
data.OPmatrix = double(data.OPall(:,8:end)); % matrix of collector readings
data.OPtime   = double(data.OPall(:,7)); % time




end % function parseTXT

