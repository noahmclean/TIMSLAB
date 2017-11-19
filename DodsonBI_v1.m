function ratiosBI = DodsonBI_v1(dataRaw, MSmethod)
%DODSONBI implements Dodson's beam interpolation algorithm as a function
%   ratiosBI is a matrix of requested beam-interpolated ratios
%   dataRaw is a matrix of non-BI intensities with columns for masses measMasses
%   MSmethod.BItimes is a vector of times in seconds for each measured mass, with settle times
%       as [peak1 settle peak2 settle peak3 settle ...]
%   MSmethod.measMasses is cell array of strings with the masses measured: {'204', '205', ...}
%   MSmethod.outRatios is a cell array of strings with the output ratios: {'204/206', '207/206', ...}
%   MSmethod.cyclesPerBlock is the number of cycles per block of data
%
%   Implementation using Ludwig's method: https://doi.org/10.1016/j.chemgeo.2009.07.004
%   from Dodson: http://iopscience.iop.org/article/10.1088/0022-3735/11/4/004/meta
%   Though note, Ludwig page 24, 2nd column should read f2 = f/2 (typo)
%
%   Noah McLean, Oct 16, 2017

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

% determine discrete times for isotope measurements in cycle
cumulativeCycleTime = cumsum(MSmethod.BItimes);
massTimes = cumulativeCycleTime(1:2:end)' - 0.5*MSmethod.BItimes(1:2:end)'; %half-way through each meas
cycleTime = cumulativeCycleTime(end);

%% Dodson setup

% Dodson BI time for each ratio (average of all four ion beam times: a1, a2, b1, b2 -> a/b)
ratioTimes = (2*massTimes(massTable(:,1)) + 2*massTimes(massTable(:,2)) + 2*cycleTime)/4;

% Ludwig's f-values - column1 for numerator, column2 for d
fValues = [ratioTimes-massTimes(massTable(:,1)) ...
                                ratioTimes-massTimes(massTable(:,2))]/cycleTime;


%% Perform beam interpolation
numBlocks = size(dataRaw,1)/cyclesPerBlock;
nBlocks = floor(numBlocks);
partialBlockCycles = max(rem(size(dataRaw,1),cyclesPerBlock)-1,0);

% create a 3D array where each row is a mass, each column a cycle, and 3rd dim = block#
dataFullBlocks = reshape(dataRaw(1:cyclesPerBlock*nBlocks,:)',[nMasses, cyclesPerBlock, nBlocks]);
firstCycle = dataFullBlocks(:,1:end-1,:);
secndCycle = dataFullBlocks(:,2:end,:);

%perform beam interpolation on complete blocks
blocksBI = ( repmat(1-fValues(:,1),1,size(firstCycle,2), size(firstCycle,3)).*firstCycle(massTable(:,1),:,:) + ...
             repmat(  fValues(:,1),1,size(secndCycle,2), size(secndCycle,3)).*secndCycle(massTable(:,1),:,:) ) ./ ...
           ( repmat(1-fValues(:,2),1,size(firstCycle,2), size(firstCycle,3)).*firstCycle(massTable(:,2),:,:) + ...
             repmat(  fValues(:,2),1,size(secndCycle,2), size(secndCycle,3)).*secndCycle(massTable(:,2),:,:) );

%perform beam interpolation on partial blocks         
if partialBlockCycles
    
    partialBlockData = dataRaw(cyclesPerBlock*nBlocks+1:end,:)';
    partialBlockBI = (repmat(1-fValues(:,1), 1, partialBlockCycles).* partialBlockData(massTable(:,1), 1:end-1) + ...
                      repmat(  fValues(:,1), 1, partialBlockCycles).* partialBlockData(massTable(:,1), 2:end) ) ./ ...
                     (repmat(1-fValues(:,2), 1, partialBlockCycles).* partialBlockData(massTable(:,2), 1:end-1) + ...
                      repmat(  fValues(:,2), 1, partialBlockCycles).* partialBlockData(massTable(:,2), 2:end) ) ;
    
else, partialBlockBI = [];
end

%% append full blocks and partial block, convert to rows = BIcycles, columns = ratios
ratiosBI = [reshape(blocksBI, [nRatios, (cyclesPerBlock-1)*nBlocks]) partialBlockBI]';

