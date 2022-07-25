function setup = sampleSetup(data)
%SAMPLESETUP set up splines, store off useful parameters
%   separate spline fit for each block
%   separate splines for 

%% spline parameters

setup.bdeg = 3; % cubic splines
setup.pord = 2; % order of differences (2nd order for min integral of 2nd derivative)
setup.scaleInt = 10; % use scaleInt as many spline coefficients as cycles
setup.scaleBeta = 1; % note different unit -- use scaleBeta betas *per block*

setup.IntLambdaInit = 1;
setup.BetaLambdaInit = 1;

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

setup.numeratorIsotopeIdx = 4; % for Sm, 150
setup.denominatorIsotopeIdx = 1; % for Sm, 147


%% constants (other than isotopic mass)

setup.kB = 1.38064852e-23; % Boltzmann's constant, J/K
setup.tempInK = 290; % temperature in K (decabin cooled to ~16 C)
setup.coulomb = 6241509074460762607.776; % 1 coulomb in elementary charges

end % function

