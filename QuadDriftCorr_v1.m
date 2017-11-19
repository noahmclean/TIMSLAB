function ratiosBI = QuadDriftCorr_v1(dataRaw, MSmethod)
%QUADDRIFCORR_v1 performs quadratic drift correction as outlined in Thermo Triton manual
%   ratiosBI is a matrix of requested beam-interpolated ratios
%   dataRaw is a matrix of non-BI intensities with columns for masses measMasses
%   MSmethod.BItimes is a vector of times in seconds for each measured mass, with settle times
%       as [peak1 settle peak2 settle peak3 settle ...]
%   MSmethod.measMasses is cell array of strings with the masses measured: {'204', '205', ...}
%   MSmethod.outRatios is a cell array of strings with the output ratios: {'204/206', '207/206', ...}
%   MSmethod.cyclesPerBlock is the number of cycles per block of data
%
%   Implemented using "In matrix notation" equation from page 2 of QuaDriftCorrection.pdf, fixing
%   typo in B(2).  Also, forcing ratio = 0 when I1(t1) = 0, which breaks intensity ratios in S and B.
%
%   Noah McLean, November 7, 2017

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

numBlocks = size(dataRaw,1)/cyclesPerBlock; %total (plus fractional) blocks of data
nBlocks = floor(numBlocks); % complete blocks of data
partialBlockCycles = max(rem(size(dataRaw,1),cyclesPerBlock)-1,0);
nCycles = nBlocks*(cyclesPerBlock-1) + partialBlockCycles;


%% Perform quadratic drift correction

S = zeros(3,3);
ratiosBI = zeros(nCycles, nRatios);

% make a vector of the first cycle of each BI pair for full blocks
cycleCount = repmat(1:(cyclesPerBlock-1),1,nBlocks);
blockTotal = repelem(0:cyclesPerBlock:(cyclesPerBlock*nBlocks-1),cyclesPerBlock - 1);
ratioStartCycles = cycleCount + blockTotal; % cycle indices for first cycle in BI ratios

% append cycles from partial blocks to end; does nothing for partialBlockCycles = 0.
ratioStartCycles = [ratioStartCycles  (1:partialBlockCycles) + cyclesPerBlock*nBlocks];

for ratioi = 1:nRatios
    
    [i1indx, minindx] = min(massTable(ratioi,:)); % i1 gets measured first
    i2indx = max(massTable(ratioi,:));            % i2 gets measured second
    
    % time between numerator and denominator isotopes
    t2 = massTimes(i2indx) - massTimes(i1indx);
    t3 = cycleTime;                               % time between first and second numerator isotope
    t4 = t2 + cycleTime;
    
    ratioCount = 1; % counter for ratios (since last cycle is skipped)
    
    for cyclei = ratioStartCycles
        
        % i1(t1) = 0 breaks ratios in S([1 3], 1) and B(2);
        if dataRaw(cyclei, i1indx) == 0
            ratioCount = ratioCount + 1;
            continue
        end 
        
        S(1,1) = dataRaw(cyclei,   i2indx) / dataRaw(cyclei, i1indx);
        S(3,1) = dataRaw(cyclei+1, i2indx) / dataRaw(cyclei, i1indx);
        S(:,2:3) = [-t2 -t2^2; t3 t3^2; -t4 -t4^2];
        B = [1
            dataRaw(cyclei+1, i1indx)/dataRaw(cyclei, i1indx) - 1
            1];
        Z = S\B;
        ratiosBI(ratioCount, ratioi) = Z(1);
        
        ratioCount = ratioCount + 1;
        
    end % for cyclei
    
    
    % if denominator isotope is measured first, use reciprocal of r.
    if minindx == 2
        nonzeroRatios = find(ratiosBI(:,ratioi)); % don't want to take reciprocal of zero values
        ratiosBI(nonzeroRatios, ratioi) = 1./ratiosBI(nonzeroRatios, ratioi);
    end % if
    
    
end % for ratioi
    
