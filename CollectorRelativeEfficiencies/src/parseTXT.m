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

% future work: change to table import for ATONA data
%opts = setvartype(opts, 'string');
%opts.DataLines = [13 Inf];
%dataTable = readtable(textFileInfo.name, opts);

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
    nDataColumns = 7 + 2*nCollectors - 1; % minus one because PM doesn't have BChannel
else
    disp('unrecognized text file column setup')
end

% grab the gains
collRange = (collectorsStartPosition+1):collectorBlockEndPosition;
data.Collectors = reshape(dtmp(collRange), 6, [])';
data.FaradayGains = double(data.Collectors(4:12,4));

% range starts after header, continues to cell before next block flag
BLrange = (baselinesStartPosition+1+nDataColumns):(onPeakStartPosition-1);
OPrange = (onPeakStartPosition+1+nDataColumns):(endPosition-1);

data.BLall = reshape(dtmp(BLrange),nDataColumns,[])';
data.OPall = reshape(dtmp(OPrange),nDataColumns,[])';

data.BLserial = double(data.BLall(:,2:4));% [block cycle integration] serially assigned counts
data.BLmatrix = double(data.BLall(:,8:end)); % matrix of collector readings
data.BLtime   = double(data.BLall(:,7)); % time
data.BLID = data.BLall(:,1); % baseline ID, eg "BL1", "BL2", etc. 1st column in TXT data file
data.BLSeqIdx = double(extractAfter(data.BLall(:,1), "BL"));

data.OPserial = double(data.OPall(:,2:4));% [block cycle integration] serially assigned counts
data.OPmatrix = double(data.OPall(:,8:end)); % matrix of collector readings
data.OPtime   = double(data.OPall(:,7)); % time
data.OPID = data.OPall(:,1); % OnPeak ID, eg "OP1", "OP2", etc.  1st column in TXT data file
data.OPSeqIdx = double(extractAfter(data.OPall(:,1), "OP"));


%% correct bad 2-second relay settle time for gains on select samples
% goodGains contains 6-second RST ccgains measured on 08-Jul-2022 (mean of 100)

if any([ "SmEfficiency_Bead3Run2-1809.TXT" "SmEfficiency_Bead3Run3-1810"] ...
               == string(methodHeader{1,2}) )
goodGains = [0.9964579, 1.0118078, 0.9885882, 1.0086657, 1.0000000, ...
             1.0149499, 0.9886156, 0.9841766, 0.9649764];

% undo bad gains, apply good gains
data.BLmatrix(:,3:11) = data.BLmatrix(:,3:11) ./ (data.FaradayGains');
data.BLmatrix(:,3:11) = data.BLmatrix(:,3:11) .* goodGains;
data.OPmatrix(:,3:11) = data.OPmatrix(:,3:11) ./ (data.FaradayGains');
data.OPmatrix(:,3:11) = data.OPmatrix(:,3:11) .* goodGains;

end % if any Sm were measured with bad 2-second RST ccgains


end % function parseTXT

