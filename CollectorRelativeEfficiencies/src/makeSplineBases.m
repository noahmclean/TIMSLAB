function B = makeSplineBases(d, setup)
%MAKESPLINEBASES Make spline bases for OP intensity and beta interpolation
%   Making these on the fly in the model evaluation took > 1 second

nBlocks = max(d.block);
nIntegPerBlock = sum(d.block == 1 & d.isOP);

B.Bint = zeros(nIntegPerBlock, setup.nCoeffInt, nBlocks);
B.Bbeta = zeros(nIntegPerBlock, setup.nCoeffBeta, nBlocks);

for iBlock = 1:nBlocks

    inBlock = d.block == iBlock & d.isOP;
    dTime = d.time(inBlock);

    Bint = bbase(dTime, ...
        setup.blockStartEndTime(iBlock,1), ...
        setup.blockStartEndTime(iBlock,2), ...
        setup.nCoeffInt-setup.bdeg, ...
        setup.bdeg);

    B.Bint(:,:,iBlock) = Bint;

    Bbeta = bbase(dTime, ...
        setup.blockStartEndTime(1,1), ...
        setup.blockStartEndTime(end,2), ...
        setup.nCoeffBeta-setup.bdeg, ...
        setup.bdeg);

    B.Bbeta(:,:,iBlock) = Bbeta;

end % for iBlock

end % function makeSplineBases