%% calculate faraday collector relative efficiencies using multidynamic measured data
% inputs: a .txt data file, and the TIMSAM method file

addpath(genpath("../data"))

%% input a filename from data folder



%% grab the corresponding methods file, make a run table for OP and BL

%method = parseTIMSAM('Sm147to150_S6.TIMSAM');
method = parseTIMSAM('Pb cup efficiency.TIMSAM');

FaraNames = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
nSeq = size(method.onpeaks, 2);
OPnames = [method.onpeaks.Name]';
nFara = size(FaraNames,2);

OPT = table('Size', [nSeq, nFara], ...
                'VariableTypes', repelem("string", nFara), ...
                'VariableNames', FaraNames, ...
                'RowNames', OPnames);
OPTable = fillmissing(OPT, 'constant', ""); 
OPMasses = fillmissing(OPT, 'constant', "0"); clear OPT

% make OP table
for iSeq = 1:nSeq

    seqName = method.onpeaks(iSeq).Name;
    seqString = string(method.onpeaks(iSeq).Info(13).Value);
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
    masses = double(extract(seqAssign(:,1), digitsPattern)');
    for iMass = 1:nMasses
        OPMasses.(activeCollectors(iMass))(seqName) = string(masses(iMass));
    end % for iMass
    
    
end % for iOP

% calculate mass differences
OPMasses = double(table2array(OPMasses));
OPMassDifTable = zeros(size(OPMasses) - [0 1]);
for iSeq = 1:nSeq
    seqMass = OPMasses(iSeq,:);
    massPresent = seqMass > 0;
    massDifPlace = massPresent(2:end) & massPresent(1:end-1);
    massVec = seqMass(massPresent);
    massDifs = massVec(2:end)-massVec(1:end-1);
    OPMassDifTable(iSeq, massDifPlace) = massDifs;
end % for iMass

OPMassDif = zeros(1, nFara-1);
for iFara = 1:nFara-1

    massDifFara = OPMassDifTable(:,iFara);
    massDifFara = massDifFara(massDifFara > 0);
    if ~all(massDifFara == massDifFara(1))
        disp('problem with OPMassDifTable')
    else 
        OPMassDif(iFara) = massDifFara(1);
    end

end % for iFara

% handle baselines
if isfield(method, 'baselines') % if baselines present

    nBL = size(method.baselines, 2);
    BLnames = [method.baselines.Name]';
    BLT = table('Size', [nBL, nFara], ...
        'VariableTypes', repelem("string", nFara), ...
        'VariableNames', FaraNames, ...
        'RowNames', BLnames);
    BLTable = fillmissing(BLT, 'constant', ""); clear BLT

end

%% 