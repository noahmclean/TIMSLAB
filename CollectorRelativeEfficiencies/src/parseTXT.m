function data = parseTXT(dataFolder, method)
%PARSETXT Parse Isotopx TXT file output
%   dataFolder is .RAW folder, method struct from parseTIMSAM

textFileInfo = dir("../data/" + dataFolder + "/*.TXT");

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

if method.settings(19).Value == "PrimaryA" % if resistor-based amplifiers, no BChannels
    nDataColumns = 7 + nCollectors;
elseif method.settings(19).Value == "SecondaryB"
    nDataColumns = 7 + 2*nCollectors;
else
    disp('unrecognized text file column setup')
end

% get method row from header, find method name
opts = delimitedTextImportOptions('NumVariables', nDataColumns);
opts.DataLines = [5, 5];
methodRow = readcell(textFileInfo.name, opts);
data.methodName = string(methodRow{2});

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

