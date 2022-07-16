%function m0 = initializeModel(data, d, method, spl)
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
%   m0.rangeRatio = indices for 'true' ratios
%   m0.rangeRefVolts = indices for reference voltages
%   m0.rangeRelEffs = indices for relative efficiencies
%   m0.rangeInts = matrix of spline coefficients for
%     intensity functions of time, each column for a block.
%   m0.rangeBetas = matrix of spline coefficients for
%     beta functions of time, each column for a block


%% establish ranges for parameters
% count up sizes: global parameters
nRatios   = 2;             % true ratios (hard-code 2 for Sm for now)
nRefVolts = 9;             % reference voltages on Faradays (BL)
nRelEffs  = nRefVolts - 1; % detector relative efficiencies

% count up sizes: blockwise parameters
nBlocks = max(data.OPserial(:,1));
nCycles = max(data.OPserial(:,2));

% scale up/down the number of spline coefficients based on cycles
scaleInt = 2; % use scaleInt as many spline coefficients as cycles
scaleBeta = 1; % use scaleBeta as many spline coeffs as cycles 

nIntCoefs  = nBlocks*nCycles*scaleInt;
nBetaCoefs = nBlocks*nCycles*scaleBeta;

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
% reshape rangeInts and Betas as a matrix, each column is a block
m0.rangeInts  = reshape(m0.rangeInts, [nCycles*scaleInt, nBlocks]);
m0.rangeBetas = reshape(m0.rangeBetas, [nCycles*scaleBeta, nBlocks]);


%% initialize global parameters

% temporary variable m will become m0.vec
m = zeros(totalModelPars,1);

% Sm ratios to 147, assume 150/147 = 1/2.031957
% Sm ratios are [148/147 149/147], 
% taken from Brennecka et al 2013 PNAS Ames Sm data (Table S5)
m(m0.rangeRatio) = [1.523370/2.031957; 1.872696/2.031957];

% reference volts from BL data, effs
meanBL = mean(data.BLmatrix)';
m(m0.rangeRefVolts) = meanBL; % ignore peak tails
m(m0.rangeRelEffs) = ones(nRelEffs,1);


%% initialize intensity functions - spline coeffs
%        L5       L4      L3     L2     Ax     H1     H2     H3     H4
%        1        2       3      4      5      6      7      8      9
%fudge = [0.0085 0.0045 -0.0045 0.0025 0.0010 0.0040 0.0003 0.0015 0.0055]';
fudge = zeros(9,1);

% set up initial i147 fit
h = axes; hold on
for monitorIsotope = [1 2 3]
%monitorIsotope = 2; % index for monitor isotope to fit
for iBlock = 7
    
    isMonitor = (d.iso == monitorIsotope) & (d.block == iBlock);
    detMonitor = d.det(isMonitor); % detector index
    intMonitor = d.int(isMonitor) - meanBL(detMonitor) + fudge(detMonitor);
    timeMonitor = d.time(isMonitor);
    [timeMonitor, sortIdx] = sort(timeMonitor);
    intMonitor = intMonitor(sortIdx);
    detMonitor = detMonitor(sortIdx);

% do some plotting for debugging
    for iDet = 1:9
        if monitorIsotope == 3, mult = 1.075; 
        elseif monitorIsotope == 2, mult =1.318;
            else, mult = 1; end
        isDet = detMonitor == iDet;
        plot(timeMonitor(isDet), mult*intMonitor(isDet), '.')

    end % for iDet

end % for iBlock

end % for monitorIsotope
%G = 


%% assign m to output structure m0

m0.vec = m;

%end % function initializeModel

