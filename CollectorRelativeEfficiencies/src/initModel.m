function m0 = initModel(data, d, tails, method, setup, B)
%INITIALIZEMODEL Initialize model vector
%   
%   model vector is (for now):
%
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

%% establish ranges for parameters
% count up sizes: global parameters
nIsotopes = max(max(method.F_ind)); % number of isotopes in method
nRatios   = nIsotopes - 2;          % number of ratios free to vary (not for norm)
nRefVolts = width(method.OPTable);  % reference voltages on Faradays (BL)
nRelEffs  = nRefVolts - 1;          % detector relative efficiencies (logged)

% count up sizes: blockwise parameters
nBlocks = max(data.OPserial(:,1));
nCycles = max(data.OPserial(:,2));

% scale up/down the number of spline coefficients based on cycles
nIntCoefs  = nBlocks*nCycles*setup.scaleInt;
nBetaCoefs = nBlocks*setup.scaleBeta;

% determine where categories start in m0
startRatio     = 1;
startRefVolt   = startRatio + nRatios;
startRelEffs   = startRefVolt + nRefVolts;
startIntCoefs  = startRelEffs + nRelEffs;
startBetaCoefs = startIntCoefs + nIntCoefs;
totalModelPars = startBetaCoefs + nBetaCoefs - 1;

% record ranges over which parameters are stored in m
m0.rangeRatio    = startRatio:(startRefVolt-1);
m0.rangeRefVolts = startRefVolt:(startRelEffs-1);
m0.rangeRelEffs  = startRelEffs:(startIntCoefs-1);
m0.rangeInts     = startIntCoefs:(startBetaCoefs-1);
m0.rangeBetas    = startBetaCoefs:totalModelPars;

% reshape rangeInts as a matrix, each column is a block
m0.rangeInts  = reshape(m0.rangeInts, [nCycles*setup.scaleInt, nBlocks]);


%% initialize global parameters

% temporary variable m will become m0.vec
m = zeros(totalModelPars,1);

% Sm ratios to 147, assume 150/147 = 1/2.031957
% Sm ratios are [148/147 149/147], 
% taken from Brennecka et al 2013 PNAS Ames Sm data (Table S5)
% all ratios are log-ratios (alr transform)
%m(m0.rangeRatio) = log([1.523370/2.031957; 1.872696/2.031957]);

freeNumeratorIdcs = 1:nIsotopes;
freeNumeratorIdcs = freeNumeratorIdcs( ...
                    freeNumeratorIdcs ~= setup.numeratorIsotopeIdx & ...
                    freeNumeratorIdcs ~= setup.denominatorIsotopeIdx);

% initialize log-ratios for "free" numerator isotopes
m(m0.rangeRatio) = log( setup.referenceMaterialIC(freeNumeratorIdcs)/...
                   setup.referenceMaterialIC(setup.denominatorIsotopeIdx) );

% assume collector relative efficiences = 1, so log(1) = 0
m(m0.rangeRelEffs) = log(ones(nRelEffs,1));


%% 1. Estimate RefVolts from BL data

% reference volts from BL data, effs
meanBL = mean(data.BLmatrix)';
m(m0.rangeRefVolts) = meanBL; % ignore peak tails

