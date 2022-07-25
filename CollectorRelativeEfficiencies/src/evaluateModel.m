function dhat = evaluateModel(d, m, m0)
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

r148147true = m(1);
r149ar7true = m(2);
refVoltages = m(m0.rangeRefVolts);
relEffs     = m(m0.rangeRelEffs);
intensities = m(m0.rangeInts(:));
intensities = reshape(intensities, size(m0.rangeInts));
betas       = m(m0.rangeBetas);

dhat = zeros(size(d.int));


%% baseline dhats

% assumption: baseline constant through analysis
BLdetectors = d.det(~d.isOP);
dhat(~d.isOP) = refVoltages(BLdetectors(BLdetectors>0));



end % evaluateModel

