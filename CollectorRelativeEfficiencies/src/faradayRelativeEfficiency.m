%% calculate faraday collector relative efficiencies using multidynamic measured data
% inputs: a .txt data file, and the TIMSAM method file

addpath(genpath("../data"))

%% 1. input a filename from data folder

%dataFolder = "Sm/SmKU1A-A2.RAW";
%dataFolder = "Sm/SmEfficiency_Bead3Run3.RAW";
%dataFolder = "Pb/21042022 NBS 982 cup efficiency.RAW"; not yet
%dataFolder = "Pb/A520_Pb.RAW";
dataFolder = "Pb/B195.RAW";

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

%% Make and store off some spline bases
% slows down evaluateModel to make them on the fly

B = makeSplineBases(d, setup);

%% 8. calculate best fit

tic
opts = optimoptions("fminunc", 'Display', 'iter-detailed');
opts.MaxFunctionEvaluations = 1e5;
[mhat, chi2] = fminunc(@(m) objfunc(d, m, s2, m0, tails, setup, B), m0.vec, opts);
toc

%% 9. estimate uncertainty in fit

dhat = evaluateModel(d, mhat, m0, tails, setup, B);
G = makeG(d, mhat, m0, tails, setup, B, dhat);
covm = inv( (G.*(1./s2))' * G ); % inv(G'*W*G), where W = diag(1./s2);
unctx = sqrt(diag(covm));
BIC = length(mhat)*log(length(dhat)) - 2*(-1/2*(chi2 + sum(log(s2))));


%% 10. visualize results

plotFit(dhat, d, mhat, m0, tails, setup, B, covm);

