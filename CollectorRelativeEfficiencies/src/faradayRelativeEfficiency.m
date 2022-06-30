%% calculate faraday collector relative efficiencies using multidynamic measured data
% inputs: a .txt data file, and the TIMSAM method file

addpath(genpath("../data"))

%% 1. input a filename from data folder

dataFolder = "Sm/SmKU1A-A2.RAW";
%dataFolder = "Pb/21042022 NBS 982 cup efficiency.RAW"; not yet

%% 2. parse the data file

data = parseTXT(dataFolder);

%% 3. grab the corresponding methods file, make run tables for OP and BL

method = parseTIMSAM(data.header.methodName);
%method = parseTIMSAM('Pb cup efficiency.TIMSAM');
%method = parseTIMSAM('PbFaraday_Pbc3Line.TIMSAM');

FaraNames = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
collectorDeltas = [-4 -3 -3 -1 0 1 2 3 4]; 
method = processMethod(method, FaraNames, collectorDeltas);
% collector mass differences (derive from OPMatrix in future)

%% 4. assemble data vector and tags

d = assembleDataVector(data, method);

%% 6. create tail model

tails = initializePeakTails(method);

%% 5. initialize model

m0 = inializeModel(data, method);


