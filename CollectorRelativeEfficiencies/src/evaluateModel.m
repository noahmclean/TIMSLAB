function dhat = evaluateModel(d, m, m0, tails, setup)
%EVALUATEMODEL Calculate dhat based on model
%   Detailed explanation goes here

%   *******From initializeModel:********
%   GLOBAL MODEL PARAMETERS:
%   148/147 true
%   149/147 true
%   L5ref reference voltage
%   ...
%   H4ref reference voltage
%   eL5   detector efficiency
%   ...
%   eH4   detector efficiency
%   
%   BLOCKWISE MODEL PARAMETERS:
%   147i(t) - intensity function for 147Sm
%   beta(t) - fractionation function
%
%   m0 is a structure that contains:
%   m0.vec - initial estimate of model vector
%   m0.rangeRatio = indices for 'true' unknown ratios
%   m0.rangeRefVolts = indices for reference voltages
%   m0.rangeRelEffs = indices for relative efficiencies
%   m0.rangeInts = matrix of spline coefficients for
%     intensity functions of time, each column for a block.
%   m0.rangeBetas = vector of spline coefficients for
%     beta function of time over the analysis


%% unpack model

nBlocks = max(d.block);

lr148147true = m(1);
lr149ar7true = m(2);
refVoltages = m(m0.rangeRefVolts);
logRelEffs  = m(m0.rangeRelEffs);
logIntsty = m(m0.rangeInts(:));
logIntsty = reshape(logIntsty, size(m0.rangeInts));
betas       = m(m0.rangeBetas);

dhat = zeros(size(d.int));

% assemble isotope ratio and isotope mass ratio vectors
% for now, hard-code Sr: [147/147, 148/147, 149/147, 150/147];
logRatioVector = [0 lr148147true lr149ar7true log(setup.internalNormRatio)];

% assemble isotope mass ratio vector needed for exponential mass
% fractionation correction
logMassVector = log([1 mass.Sm148/mass.Sm147 mass.Sm149/mass.Sm147 mass.Sm150/mass.Sm147]);


%% baseline dhats

% assumption: reference voltage constant through analysis
% assumption: peak tails for BL are relative to first intensity measured
%             during the block (taken as first knot of each block for
%             spline fit to major peak intensity)

isBL = ~d.isOP;

% get detector vector, just for BL integrations
BLdetectors = d.det(isBL); 

% reference voltages, assigned by detector, for each integration
BLrefVolts  = refVoltages(BLdetectors(BLdetectors>0)); 

% peak tails for BL integrations
BLblocks = d.block(isBL);
majorPeakIntAtBlockStarts = exp(logIntsty(1,:));
BLmajorPeakInt = majorPeakIntAtBlockStarts(BLblocks)';
% peak tails scaled to major isotope intensity:
BLpeakTailSums = tails.dvec(isBL) .* BLmajorPeakInt; 

% baseline meas = ref voltage + peak tail
dhat(isBL) = BLrefVolts + BLpeakTailSums;


%% on-peak dhats: reference voltages

isOP = d.isOP;

% get detector vector, just for OP integrations
OPdetectors = d.det(isOP); 

% reference voltages, assigned by detector, for each integration
OPrefVolts  = refVoltages(OPdetectors(OPdetectors>0)); 


%% on-peak dhats: fit intensities

for iBlock = 1:nBlocks

    inBlock = d.block == iBlock & d.isOP;
    dTime = d.time(inBlock);
    dInt  = d.int(inBlock);
    dIso = d.iso(inBlock);

    % primary beam intensity
    Bint = bbase(dTime, ...
              setup.blockStartEndTime(iBlock,1), ... 
              setup.blockStartEndTime(iBlock,2), ...
              setup.nCoeffInt-setup.bdeg, ...
              setup.bdeg);

    knotsForiBlock = logIntsty(:,iBlock);
    primaryBeamLogInt = Bint*knotsForiBlock;

    % fit for beta. note: (conventaional beta)/denominatorMass, units /amu
    % 
    Bbeta = bbase(dTime, ...
              setup.blockStartEndTime(1,1), ... 
              setup.blockStartEndTime(end,2), ...
              setup.nCoeffBeta-setup.bdeg, ...
              setup.bdeg);
    betaForiBlock = Bbeta*betas;


    

    

end % for iBLock




end % evaluateModel

