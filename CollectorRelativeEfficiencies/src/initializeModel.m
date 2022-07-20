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
%   m0.rangeRatio = indices for 'true' unknown ratios
%   m0.rangeRefVolts = indices for reference voltages
%   m0.rangeRelEffs = indices for relative efficiencies
%   m0.rangeInts = matrix of spline coefficients for
%     intensity functions of time, each column for a block.
%   m0.rangeBetas = matrix of spline coefficients for
%     beta function of time, each column for a block


%% establish ranges for parameters
% count up sizes: global parameters
nRatios   = 2;             % true ratios (hard-code 2 for Sm for now)
nRefVolts = 9;             % reference voltages on Faradays (BL)
nRelEffs  = nRefVolts - 1; % detector relative efficiencies

% count up sizes: blockwise parameters
nBlocks = max(data.OPserial(:,1));
nCycles = max(data.OPserial(:,2));

% scale up/down the number of spline coefficients based on cycles
nIntCoefs  = nBlocks*nCycles*spl.scaleInt;
nBetaCoefs = nBlocks*nCycles*spl.scaleBeta;

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
m0.rangeInts  = reshape(m0.rangeInts, [nCycles*spl.scaleInt, nBlocks]);
m0.rangeBetas = reshape(m0.rangeBetas, [nCycles*spl.scaleBeta, nBlocks]);


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

% set up initial i147 fit, with nSegInt knots
monitorIsotope = 2; % index for monitor isotope to fit
D = diff(eye(spl.nSegInt), spl.pord); % 2nd order smoothing, cubic spline;
for iBlock = 1:nBlocks
    
    isMonitor = (d.iso == monitorIsotope) & (d.block == iBlock);
    detMonitor = d.det(isMonitor); % detector index
    intMonitor = d.int(isMonitor) - meanBL(detMonitor);
    timeMonitor = d.time(isMonitor);
    [timeMonitor, sortIdx] = sort(timeMonitor);
    intMonitor = intMonitor(sortIdx);
    detMonitor = detMonitor(sortIdx);

    % set up spline basis
    B = bbase(timeMonitor, min(timeMonitor), max(timeMonitor), ...
              spl.nSegInt-spl.bdeg, spl.bdeg);
    lambda = spl.IntLambdaInit;
    Baugmented = [B; sqrt(lambda)*D];
    yaugmented = [intMonitor; zeros(size(D,1),1)];
    
    % least squares fit, no weights
    %xls = (B'*B)\(B'*intMonitor);
    %plot(timeMonitor, intMonitor, '.'); hold on
    %plot(timeMonitor, B*xls, '-k')

    % smoothing spline
    xsmooth = (Baugmented'*Baugmented)\(Baugmented'*yaugmented);
    %plot(timeMonitor, intMonitor, '.'); hold on
    %plot(timeMonitor, B*xsmooth, '-r')

    m(m0.rangeInts(:,iBlock)) = xsmooth;

end % for iBlock


%% initialize spline coefficients -- betas

cmap = lines(9);
cmap(8,:) = [1 1 1];
cmap(9,:) = [1 0 0];
numeratorIsotopeIdx = 4;
denominatorIsotopeIdx = 1;

isNumerator = (d.iso == numeratorIsotopeIdx);
isDenominator = (d.iso == denominatorIsotopeIdx);
timeNumerator = d.time(isNumerator);
timeDenominator = d.time(isDenominator);
isMeasBoth = timeNumerator == timeDenominator;
detNumer = d.det(isNumerator);
detDenom = d.det(isDenominator);
numIntensity = d.int(isNumerator) - meanBL(detNumer);
denIntensity = d.int(isDenominator) - meanBL(detDenom);
numIntensity = numIntensity(isMeasBoth);
denIntensity = denIntensity(isMeasBoth);
timeRatio = timeNumerator(isMeasBoth);





%% assign m to output structure m0

m0.vec = m;

%end % function initializeModel

