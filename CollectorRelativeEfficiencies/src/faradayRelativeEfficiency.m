%% calculate faraday collector relative efficiencies using multidynamic measured data
% inputs: a .txt data file, and the TIMSAM method file

addpath(genpath("../data"))

%% input a filename from data folder



%% grab the corresponding methods file, make a run table for OP and BL

method = parseTIMSAM('Sm147to150_S6.TIMSAM');

FaraNames = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
nBL = size(method.baselines, 2);
BLnames = [method.baselines.Name]';
nOP = size(method.onpeaks, 2);
OPnames = [method.onpeaks.Name]';
nFara = size(FaraNames,2);

BLTable = table('Size', [nBL, nFara], ...
                'VariableTypes', repelem("string", nFara), ...
                'VariableNames', FaraNames, ...
                'RowNames', BLnames);
BLTable = fillmissing(BLTable, 'constant', "");
OPTable = table('Size', [nOP, nFara], ...
                'VariableTypes', repelem("string", nFara), ...
                'VariableNames', FaraNames, ...
                'RowNames', OPnames);
OPTable = fillmissing(OPTable, 'constant', "");

% make OP table
for iOP = 1:nOP

    seqName = method.onpeaks(iOP).Name;
    seqString = string(method.onpeaks(iOP).Info(13).Value);
    seqAssign = split(seqString, ",");
    seqAssign = split(seqAssign, ":");

    % detectors names are concatenated with sequence name, eg "H3S1"
    activeCollectors = extractBefore(seqAssign(:,2), seqName)';
    
    % assign massID to OP sequence table
    nMasses = size(seqAssign, 1);
    for iMass = 1:nMasses
        OPTable.(activeCollectors(iMass))(seqName) = seqAssign(iMass,1);
    end

    % record approx isotopic masses (N values) in mass table for deltas
    

end % for iOP

% calculate mass differences



%% 