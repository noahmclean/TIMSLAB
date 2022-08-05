function [] = plotFit(dhat, d, mhat, m0, tails, setup, B, covx)
%PLOTFIT Plots to assess model fit
%   Plot 1: Baselines
%   Plot 2: Beam interpolation on major peak
%   Noah McLean for Faraday Relative Efficiency project, 31 July 2022

reduceData = @(reductionFactor, dataset) ...
              transpose(sum(reshape(dataset, reductionFactor, []), 1) / reductionFactor);

%% Plot 1: Baselines

f1 = figure('Position', [1 1 2000 1200], 'Units', 'pixels');
isBL = ~d.isOP;
refVoltsAll = mhat(m0.rangeRefVolts);

for iDet = 1:9

    subplot(3,3,iDet); hold on
    BLfromDet = isBL & d.det == iDet;
    BLmeas = d.int(BLfromDet);
    tailsForDet = tails.dvec(BLfromDet);
    timeForDet = d.time(BLfromDet);
    [~, idxSort] = sort(timeForDet);
    
    refVoltsFit = refVoltsAll(iDet);

    % reduce data size for plotting by factor reductionFactorForPlottingBL

    reductionFactor = setup.reductionFactorForPlottingBL;
    
    idxSortReduced = reduceData(reductionFactor, idxSort);
    BLmeasReduced = reduceData(reductionFactor, BLmeas(idxSort));
    tailsReduced = reduceData(reductionFactor, tailsForDet(idxSort));

    plot(idxSortReduced, BLmeasReduced, '.b')
    plot(idxSortReduced, refVoltsFit, '.r')
    plot(idxSortReduced, refVoltsFit + tailsReduced, '-g')


end % for iDet, baselines

set(f1, 'Name', 'Baselines and Peak Tails')


%% Plot 2: Major isotope beam fit

nBlocks = max(d.block);

lr148147true = mhat(1);
lr149ar7true = mhat(2);
refVoltages = mhat(m0.rangeRefVolts);

logRelEffs = mhat(m0.rangeRelEffs);
logRelEffs = [logRelEffs(1:4); 0; logRelEffs(5:8)]; % add axial... fix later

logIntensity = mhat(m0.rangeInts(:));
logIntensity = reshape(logIntensity, size(m0.rangeInts));
betas       = mhat(m0.rangeBetas); % in per amu units

logRatioVector = [0 lr148147true lr149ar7true log(setup.internalNormRatio)]';

logMassRatioVector = ...
              log([1 mass.Sm148/mass.Sm147 mass.Sm149/mass.Sm147 mass.Sm150/mass.Sm147])';

denomMass = mass.Sm147;

for iBlock = 1:nBlocks

    inBlock = d.block == iBlock & d.isOP;
    dDet  = d.det(inBlock);
    dIso = d.iso(inBlock);
    dInt = d.int(inBlock);
    dTime = d.time(inBlock);

    knotsForiBlock = logIntensity(:,iBlock);
    beamUnique = B.BintUnique(:,:,iBlock)*knotsForiBlock;
    primaryBeamLogInt = beamUnique(B.uniqueIdx(:,iBlock));


    betaUnique = B.BbetaUnique(:,:,iBlock)*betas;
    betaForiBlock = betaUnique(B.uniqueIdx(:,iBlock));

    logIsotopeRatio = logRatioVector(dIso);

    logMassRatio = logMassRatioVector(dIso);

    refVolts = refVoltages(dDet);

    tailContribution = exp(primaryBeamLogInt) .* tails.dvec(inBlock);

    logEfficiencyRatio = logRelEffs(dDet);

    corrMeas = (dInt - refVolts) ./ exp(logEfficiencyRatio);
    corrFit  = exp(logIsotopeRatio + betaForiBlock .* logMassRatio*denomMass + ...
        primaryBeamLogInt) + tailContribution;

    figure()
    plot(dTime, corrMeas, '.b'); hold on
    plot(dTime, corrFit, '.r')

end % for iBlock (plot 2)


end % function plotFit



