function [ratiosBI, blockstats] = noBI(dataRaw, MSmethod)
%NOBI calculates ratios for data that has already been beam-interpolated
%
%   dataRaw is a matrix of non-BI intensities with columns for masses measMasses
%
%   MSmethod is a structure with the following mass spectrometer method information:
%   MSmethod.BItimes is a vector of times in seconds for each measured mass, with settle times
%       as [peak1 settle peak2 settle peak3 settle ...]
%   MSmethod.measMasses is cell array of strings with the masses measured: {'204', '205', ...}
%   MSmethod.outRatios is a cell array of strings with the output ratios: {'204/206', '207/206', ...}
%   MSmethod.cyclesPerBlock is the number of cycles per block of data
%
%   Noah McLean, May 30, 2018

%% Check inputs

narginchk(2,2) % check for both arguments

cyclesPerBlock = MSmethod.cyclesPerBlock;

if ~isa(dataRaw,'double')
  error('Error: Raw data must be numeric matrix of ion beam intensities')
end

if ~isa(MSmethod.BItimes,'double') || min(size(MSmethod.BItimes)) ~= 1
  error('Error: BItimes must be numeric vector of ion beam integration times')
end

if length(MSmethod.BItimes) ~= size(dataRaw,2)*2
  error('Error: BItimes must have an integration and settle time for each isotope in dataRaw')
end

if ~isa(MSmethod.measMasses,'cell')
  error('Error: measMasses should be a cell array of strings')
end

if ~isa(MSmethod.outRatios,'cell')
  error('Error: outRatios should be a cell array of strings')
end


%% Parse info on inputs and outputs

% determine which masses go in which ratios: create a table massTable with a row for each 
% ratio and column1 = index(numerator), column2 = index(denominator)

nMasses = size(dataRaw,2);
nRatios = size(MSmethod.outRatios,2);

massInRatioMatrix = zeros(nRatios, nMasses);
massTable = zeros(nRatios,2);
for ii = 1:nRatios
    for jj = 1:nMasses
        
        massInRatioPosition = strfind(MSmethod.outRatios{ii}, MSmethod.measMasses{jj});
        if ~isempty(massInRatioPosition)
            massInRatioMatrix(ii,jj) = massInRatioPosition;
        end
        
    end
    
    massTable(ii, 1) = find(massInRatioMatrix(ii,:) == 1);
    massTable(ii, 2) = find(massInRatioMatrix(ii,:) > 1);
    
end


%% Calculate ratios

ratiosBI = dataRaw(:,massTable(:,1))./dataRaw(:,massTable(:,2));