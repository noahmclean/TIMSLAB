function setup = sampleSetup(data, method)
%SAMPLESETUP set up splines and run parametesr, store off useful parameters
%   separate spline fit for each block
%   one spline fit across all betas
%   use method to figure out which parameters to use (generalize later?)

%% spline parameters

switch method.methodName
    case "Sm147to150_S6_v2" || "Sm147to150_S6"

        setup.bdeg = 3; % 1 for linear, 2 for quadratic, 3 for cubic splines
        setup.pord = 2; % order of differences (2nd order for min integral of 2nd derivative)
        setup.scaleInt = 10; % use scaleInt as many spline coefficients as cycles
        setup.scaleBeta = 1; % note different unit -- use scaleBeta betas *per block*

        setup.IntLambdaInit = 1e-8;
        setup.BetaLambdaInit = 1;

    case "Pb_Faraday_MD_981"

        setup.bdeg = 3; % 1 for linear, 2 for quadratic, 3 for cubic splines
        setup.pord = 2; % order of differences (2nd order for min integral of 2nd derivative)
        setup.scaleInt = 1; % use scaleInt as many spline coefficients as cycles
        setup.scaleBeta = 1; % note different unit -- use scaleBeta betas *per block*

        setup.IntLambdaInit = 1e-8;
        setup.BetaLambdaInit = 1;

end

%% block start/stop and parameter ranges

nBlocks = max(data.OPserial(:,1));
nCycles = max(data.OPserial(:,2));

setup.nCoeffInt  = ceil(nCycles*setup.scaleInt);
setup.nCoeffBeta = ceil(nBlocks*setup.scaleBeta);

% block start and stop indices, times
blockStartEndIdx  = zeros(nBlocks,2);
blockStartEndTime = zeros(nBlocks,2);
for iBlock = 1:nBlocks

    blockStartEndIdx(iBlock,1) = find(data.OPserial(:,1) == iBlock, 1, 'first');
    blockStartEndIdx(iBlock,2) = find(data.OPserial(:,1) == iBlock, 1, 'last');
    blockStartEndTime(iBlock,:)= data.OPtime(blockStartEndIdx(iBlock,:))';

end % for iBlock

setup.blockStartEndIdx = blockStartEndIdx;
setup.blockStartEndTime = blockStartEndTime;


%% internal normalization
% assume denominator isotope is major isotope

switch method.methodName
    case "Sm147to150_S6_v2" || "Sm147to150_S6"

        setup.numeratorIsotopeIdx = 4; % for Sm, 150
        setup.denominatorIsotopeIdx = 1; % for Sm, 147
        % taken from Brennecka et al 2013 PNAS Ames Sm data (Table S5)
        % internalNormRatio = numeratorIsotopeIdx/denominatorIsotopeIdx
        setup.internalNormRatio = 1/2.031957;

    case "Pb_Faraday_MD_981"

        setup.numeratorIsotopeIdx = 2; % for Pb, 206
        setup.denominatorIsotopeIdx = 4; % for Pb, 208
        % taken from Brennecka et al 2013 PNAS Ames Sm data (Table S5)
        % internalNormRatio = numeratorIsotopeIdx/denominatorIsotopeIdx
        setup.internalNormRatio = 2.1681;

end % switch case methodName

%% first finite difference options

setup.relativeStepSize = 1e-7;


%% plotting options
setup.reductionFactorForPlottingBL = 100;


%% constants (other than isotopic mass)

setup.kB = 1.38064852e-23; % Boltzmann's constant, J/K
setup.tempInK = 290; % temperature in K (decabin cooled to ~16 C)
setup.coulomb = 6241509074460762607.776; % 1 coulomb in elementary charges

setup.noiseConstantATONAsCPS = 1484; % cps noise
setup.noiseATONAs = 

end % function

