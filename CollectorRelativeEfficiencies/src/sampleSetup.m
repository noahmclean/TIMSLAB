function setup = sampleSetup(data, method)
%SAMPLESETUP set up splines and run parametesr, store off useful parameters
%   separate spline fit for each block
%   one spline fit across all betas
%   use method to figure out which parameters to use (generalize later?)

%% spline parameters, depend on method (not on reference material)

if any(method.methodName == [ "Sm147to150_S6_v2" ,"Sm147to150_S6"])

    setup.bdeg = 3; % 1 for linear, 2 for quadratic, 3 for cubic splines
    setup.pord = 2; % order of differences (2nd order for min integral of 2nd derivative)
    setup.scaleInt = 10; % use scaleInt as many spline coefficients as cycles
    setup.scaleBeta = 1; % note different unit -- use scaleBeta betas *per block*

    setup.IntLambdaInit = 1e-8;
    setup.BetaLambdaInit = 1;

elseif method.methodName == "Pb_Faraday_MD_981"

    setup.bdeg = 3; % 1 for linear, 2 for quadratic, 3 for cubic splines
    setup.pord = 2; % order of differences (2nd order for min integral of 2nd derivative)
    setup.scaleInt = 1; % use scaleInt as many spline coefficients as cycles
    setup.scaleBeta = 1; % note different unit -- use scaleBeta betas *per block*

    setup.IntLambdaInit = 1e-8;
    setup.BetaLambdaInit = 1;

elseif method.methodName == "PbFaraday_Pbc3Line"

    setup.bdeg = 3; % 1 for linear, 2 for quadratic, 3 for cubic splines
    setup.pord = 2; % order of differences (2nd order for min integral of 2nd derivative)
    setup.scaleInt = 3; % use scaleInt as many spline coefficients as cycles
    setup.scaleBeta = 1; % note different unit -- use scaleBeta betas *per block*

    setup.IntLambdaInit = 1e-8;
    setup.BetaLambdaInit = 1;

elseif method.methodName == "PbFaraday_981_6Seq"

    setup.bdeg = 3; % 1 for linear, 2 for quadratic, 3 for cubic splines
    setup.pord = 2; % order of differences (2nd order for min integral of 2nd derivative)
    setup.scaleInt = 10; % use scaleInt as many spline coefficients as cycles
    setup.scaleBeta = 1; % note different unit -- use scaleBeta betas *per block*

    setup.IntLambdaInit = 1e-18;
    setup.BetaLambdaInit = 1;

else % method name not yet set up

    disp("Method not yet set up in sampleSetup.m")
    return

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


%% internal normalization and reference material IC
% assume denominator isotope is major isotope

if any(method.methodName == [ "Sm147to150_S6_v2" ,"Sm147to150_S6"])

        setup.numeratorIsotopeIdx = 4; % for Sm, 150
        setup.denominatorIsotopeIdx = 1; % for Sm, 147, major isotope of 147-150
        setup.numeratorMassID = "150Sm";
        setup.denominatorMassID = "147Sm";
        % taken from Brennecka et al 2013 PNAS Ames Sm data (Table S5)
        % internalNormRatio = numeratorIsotopeIdx/denominatorIsotopeIdx
        setup.referenceMaterialIC = [2.031957 1.523370 1.872696 1];


elseif method.methodName == "Pb_Faraday_MD_981"  || ... or
       method.methodName == "PbFaraday_Pbc3Line" || ... or
       method.methodName == "PbFaraday_981_6Seq"

        setup.numeratorIsotopeIdx = 2; % for Pb, 206
        setup.denominatorIsotopeIdx = 4; % for Pb, 208, major isotope
        setup.numeratorMassID = "206Pb";
        setup.denominatorMassID = "208Pb";
        setup.referenceMaterialIC = [0.0590074 1 0.914683 2.1681]; % Condon

end % switch case methodName

setup.internalNormRatio = setup.referenceMaterialIC(setup.numeratorIsotopeIdx)/...
                                  setup.referenceMaterialIC(setup.denominatorIsotopeIdx);

% MassName = MassID with letters and numbers switched, to index into mass class
setup.numeratorMassName   = extract(setup.numeratorMassID,  lettersPattern) + ...
                            extract(setup.numeratorMassID,  digitsPattern); 
setup.denominatorMassName = extract(setup.denominatorMassID,lettersPattern) + ...
                            extract(setup.denominatorMassID,digitsPattern); 


%% first finite difference options

setup.relativeStepSize = 1e-7;


%% plotting options

if any(method.methodName == [ "Sm147to150_S6_v2" ,"Sm147to150_S6"])

    setup.reductionFactorForPlottingBL = 100;
    
elseif method.methodName == "Pb_Faraday_MD_981" || ...
       method.methodName == "PbFaraday_Pbc3Line"

    setup.reductionFactorForPlottingBL = 8;

elseif method.methodName == "PbFaraday_981_6Seq"
    
    setup.reductionFactorForPlottingBL = 1;

end % switch case methodName



%% constants (other than isotopic mass)

setup.kB = 1.38064852e-23; % Boltzmann's constant, J/K
setup.tempInK = 290; % temperature in K (decabin cooled to ~16 C)
setup.coulomb = 6241509074460762607.776; % 1 coulomb in elementary charges

% Ickert: "1/t * 1484 = amplifier noise, t in seconds, noise in cps"
setup.noiseConstantATONAsCPS = 1484; % cps noise
R = 10^11; % software reports volts consistent with 10^11 ohm resistor
cpsPerVolt = setup.coulomb/R; % 1 coulomb in 'atomic units'.
% variance for atonas is 1/t * setup.noiseATONAs, where t is in seconds.
setup.noiseATONAs = setup.noiseConstantATONAsCPS / cpsPerVolt;


end % function

