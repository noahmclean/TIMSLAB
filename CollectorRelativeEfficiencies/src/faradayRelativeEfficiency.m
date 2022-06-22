%% calculate faraday collector relative efficiencies using multidynamic measured data
% inputs: a .txt data file, and the TIMSAM method file

addpath(genpath("../data"))

%% input a filename from data folder

dataFolder = "Sm/SmKU1A-A2.RAW";


%% grab the corresponding methods file, make run tables for OP and BL

method = parseTIMSAM('Sm147to150_S6.TIMSAM');
%method = parseTIMSAM('Pb cup efficiency.TIMSAM');
%method = parseTIMSAM('PbFaraday_Pbc3Line.TIMSAM');

FaraNames = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
method = processMethod(method, FaraNames);


%% parse the data file

data = parseTXT(dataFolder, method);
