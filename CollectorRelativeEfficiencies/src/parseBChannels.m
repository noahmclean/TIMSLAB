function [data, dataB] = parseBChannels(data)
%PARSEBCHANNELS Condense ATONA data from B Channels 
%   New channel dataB is same shape as data, but with BChannel data.
%   Bchannel data deleted from data struct, now only has AChannel data.

dataB = data; % make a copy

% remove missing values from BL data matrices
isBLData = ~isnan(data.BLmatrixB(:,1));
dataB.BLmatrix = data.BLmatrixB(isBLData,:);
dataB.BLserial = data.BLserial(isBLData,:);
dataB.BLall    = data.BLall(isBLData,:);

% future: use average cycle time instead of last time. 
% useful if cycle times differ among sequences
dataB.BLtime = data.BLtime(isBLData,:);

dataB.BLID = data.BLID(isBLData,:);
dataB.BLSeqIdx = data.BLSeqIdx(isBLData,:);

% remove missing values from OP data matrices
isOPData = ~isnan(data.OPmatrixB(:,1));

dataB.OPmatrix = data.OPmatrixB(isOPData,:);
dataB.OPserial = data.OPserial(isOPData,:);
dataB.OPall    = data.OPall(isOPData,:);

% future: use average cycle time instead of last time. 
% useful if cycle times differ among sequences
dataB.OPtime = data.OPtime(isOPData,:);

dataB.OPID = data.OPID(isOPData,:);
dataB.OPSeqIdx = data.OPSeqIdx(isOPData,:);

% assign A channel data matrices to data
data.BLmatrix = data.BLmatrixA;
data.OPmatrix = data.OPmatrixA;

% cleanup in aisle B
fields = {'OPmatrixA', 'OPmatrixB', 'BLmatrixA', 'BLmatrixB'};
data = rmfield(data,fields);
dataB = rmfield(dataB, fields);

end

