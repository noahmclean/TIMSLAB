function m0 = initializeModel(data, d, method, setup, Bstruct)
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

% reference volts from BL data, effs
meanBL = mean(data.BLmatrix)';
m(m0.rangeRefVolts) = meanBL; % ignore peak tails

% assume collector relative efficiences = 1, so log(1) = 0
m(m0.rangeRelEffs) = log(ones(nRelEffs,1));


%% initialize intensity functions - spline coeffs

% set up initial i147 fit, with nSegInt knots
% note: fit is to log-intensity of major/monitor isotope
monitorIsotope = setup.denominatorIsotopeIdx; % index for monitor isotope to fit
DInt = diff(eye(setup.nCoeffInt), setup.pord); % 2nd order smoothing, cubic spline;
for iBlock = 1:nBlocks
    
    isMonitor = (d.iso == monitorIsotope) & (d.block == iBlock);
    detMonitor = d.det(isMonitor); % detector index
    intMonitor = d.int(isMonitor) - meanBL(detMonitor); % rough BL correction
    timeMonitor = d.time(isMonitor);
    
    % sort by time (originally sorted by detector)
    [timeMonitor, sortIdx] = sort(timeMonitor);
    intMonitor = intMonitor(sortIdx);
    intMonitor = log(intMonitor); % fit (true) log-intensities, always
    % detMonitor = detMonitor(sortIdx); % may need at some point?

    % set up spline basis
    B = bbase(timeMonitor, ...
              setup.blockStartEndTime(iBlock,1), ... 
              setup.blockStartEndTime(iBlock,2), ...
              setup.nCoeffInt-setup.bdeg, ...
              setup.bdeg);
    lambda = setup.IntLambdaInit;
    Baugmented = [B; sqrt(lambda)*DInt];
    yaugmented = [intMonitor; zeros(size(DInt,1),1)];
    
    % NEW
    B = Bstruct.BintUnique(:,:,iBlock);
    lambda = setup.IntLambdaInit;
    Baugmented = [B; sqrt(lambda)*DInt];
    yaugmented = [intMonitor; zeros(size(DInt,1),1)];

    % least squares fit, no weights
%     xls = (B'*B)\(B'*intMonitor);
%     plot(timeMonitor, intMonitor, '.'); hold on
%     plot(timeMonitor, B*xls, '-k')

    % smoothing spline
    ysmooth = (Baugmented'*Baugmented)\(Baugmented'*yaugmented);
     plot(timeMonitor, exp(intMonitor), '.b', 'MarkerSize', 15); %hold on
     plot(timeMonitor, exp(B*ysmooth), '-+r', 'LineWidth', 2)

    m(m0.rangeInts(:,iBlock)) = ysmooth;    


end % for iBlock


%% initialize spline coefficients -- betas
% fit one spline across all blocks of analysis
% note: current code assumes num and denom always appear in same
% sequences

% pick out numerator and denominator, correct for baselines
isNumerator = (d.iso == setup.numeratorIsotopeIdx);
isDenominator = (d.iso == setup.denominatorIsotopeIdx);
detNumer = d.det(isNumerator);
detDenom = d.det(isDenominator);
numIntensity = d.int(isNumerator) - meanBL(detNumer);
denIntensity = d.int(isDenominator) - meanBL(detDenom);

% make sure times for numerator and denominator match
% check this for squences without both isotopes?
% timeNumerator = d.time(isNumerator);
% timeDenominator = d.time(isDenominator);
% isMeasBoth = timeNumerator == timeDenominator;
% numIntensity = numIntensity(isMeasBoth);
% denIntensity = denIntensity(isMeasBoth);
% timeRatio = timeNumerator(isMeasBoth);

isMeasBoth = isNumerator & isDenominator;
numIntensity = numIntensity(isMeasBoth);
denIntensity = denIntensity(isMeasBoth);
timeRatio = d.time(isMeasBoth);

% sort by time (originally sorted by detector)
[timeRatio, sortIdx] = sort(timeRatio);
numIntensity = numIntensity(sortIdx);
denIntensity = denIntensity(sortIdx);

% set up spline basis
BBeta = bbase(timeRatio, ...
              setup.blockStartEndTime(1,1), ...
              setup.blockStartEndTime(end,2), ...
              setup.nCoeffBeta-setup.bdeg, ...
              setup.bdeg);
lambda = setup.BetaLambdaInit;
DBeta = diff(eye(setup.nCoeffBeta), setup.pord); % 2nd order smoothing, cubic spline;
Baugmented = [BBeta; sqrt(lambda)*DBeta];

% calculate betas, using "fractionation per amu" version of formula
% beta/m147 = (log(150/147)meas - log(150/147)true)/(log(M150/M147)*m147)
% call the beta/m147 beta for simplicity

beta = (log(numIntensity./denIntensity) - log(setup.internalNormRatio)) / ...
       (log(mass.(setup.numeratorMassName)/mass.(setup.denominatorMassName)) * ...
                                                 mass.(setup.denominatorMassName));

yaugmented = [beta; zeros(size(DBeta,1),1)];
ysmooth = (Baugmented'*Baugmented)\(Baugmented'*yaugmented);

figure; hax = axes(); hold on
plot(timeRatio, beta, '.', 'MarkerSize', 5)
plot(timeRatio, BBeta*ysmooth, '-r')

m(m0.rangeBetas) = ysmooth;


%% assign m to output structure m0

m0.vec = m;


end % function initializeModel

