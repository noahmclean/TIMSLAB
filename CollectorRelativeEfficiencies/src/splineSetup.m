function spl = splineSetup(data)
%SPLINESETUP set up splines, store off useful parameters
%   separate spline fit for each block
%   separate splines for 

spl.bdeg = 3; % cubic splines
scaleInt = 2; % use scaleInt as many spline coefficients as cycles
scaleBeta = 1; % use scaleBeta as many spline coeffs as cycles 

%% setup

nBlocks = max(data.OPserial(:,1));
nCycles = max(data.OPserial(:,2));

spl.nSegInt  = nCycles*scaleInt;
spl.nSegBeta = nCycles*scaleBeta;

% block start and stop indices, times
blockStartEndIdx  = zeros(2,nBlocks);
blockStartEndTime = zeros(2,nBlocks);
for iBlock = 1:nBlocks

    blockStartEndIdx(1,iBlock) = find(data.OPserial(:,1) == iBlock, 1, 'first');
    blockStartEndIdx(2,iBlock) = find(data.OPserial(:,1) == iBlock, 1, 'last');
    blockStartEndTime(:,iBlock)= data.OPtime(blockStartEndIdx(:,iBlock));

end % for iBlock

spl.blockStartEndIdx = blockStartEndIdx;
spl.blockStartEndTime = blockStartEndTime;

end % function

