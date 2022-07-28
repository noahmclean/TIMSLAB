%% calculate faraday collector relative efficiencies using multidynamic measured data
% inputs: a .txt data file, and the TIMSAM method file

addpath(genpath("../data"))

%% 1. input a filename from data folder

%dataFolder = "Sm/SmKU1A-A2.RAW";
dataFolder = "Sm/SmEfficiency_Bead3Run3.RAW";
%dataFolder = "Pb/21042022 NBS 982 cup efficiency.RAW"; not yet
%dataFolder = "Pb/A520_Pb.RAW";

%% 2. parse the data file

data = parseTXT(dataFolder);

%% 3. grab the corresponding methods file, make run tables for OP and BL

method = parseTIMSAM(data.header.methodName);
%method = parseTIMSAM('Pb cup efficiency.TIMSAM');
%method = parseTIMSAM('PbFaraday_Pbc3Line.TIMSAM');

FaraNames = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
collectorDeltas = [-4 -3 -2 -1 0 1 2 3 4]; 
method = processMethod(method, FaraNames, collectorDeltas);
% collector mass differences (derive from OPMatrix in future)



%% 4. assemble data vector and tags

[d, data] = assembleDataVector(data, method);
% also trims data matrices to collectors from method

%% 5. create tail model

tails = initializePeakTails(method, d);

%% 5b. create spline setup

% contains nseg, block start/stop times, nSeg
setup = sampleSetup(data);

%% 6. initialize model

m0 = initializeModel(data, d, setup);

%% 7. calculate uncertainty in data

s2 = calculateUnct(data, d, method, setup);

%% timing

tic
chi2 = objfunc(d, m0.vec, s2, m0, tails, setup);
toc


%% 8. calculate best fit

opts = optimoptions("fminunc", 'Display', 'iter-detailed');
opts.StepTolerance = 1e-6;
[mhat, chi2] = fminunc(@(m) objfunc(d, m, s2, m0, tails, setup), m0.vec, opts);

% %% 9. estimate uncertainty in fit
% 
% G = makeG(d, mhat);
% covx = inv( (G.*(1./s2))' * G ); % inv(G'*W*G), where W = diag(1./s2);
% unctx = sqrt(diag(covx));
% 
% %% 10. visualize results