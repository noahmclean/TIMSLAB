function method = processMethod(method, FaraNames)
%PROCESSMETHOD Extract useful information from parsed method structure
%
%   Then append the useful info to the method structure

nSeq = size(method.onpeaks, 2);
OPnames = [method.onpeaks.Name]';
nFara = size(FaraNames,2);

% determine collector in axial position (Axial Faraday or Daly/Photomultiplier)
if string(method.settings(11).Value) == "Axial"
    method.axialPositionDetector.Name = "Axial";
    method.axialPositionDetector.Code = "Ax";
elseif string(method.settings(11).Value) == "PhotoMultiplier"
    method.axialPositionDetector.Name = "PhotoMultiplier";
    method.axialPositionDetector.Code = "PM";
else
    disp("Axial Position Detector Not Recognized!!")
end

OPT = table('Size', [nSeq, nFara], ...
                'VariableTypes', repelem("string", nFara), ...
                'VariableNames', FaraNames, ...
                'RowNames', OPnames);
OPTable = fillmissing(OPT, 'constant', ""); 
OPMasses = fillmissing(OPT, 'constant', "NaN"); clear OPT

% make OP table
for iSeq = 1:nSeq

    seqName = method.onpeaks(iSeq).Name;
    seqString = string(method.onpeaks(iSeq).Info(13).Value);
    seqAssign = split(seqString, ",");
    seqAssign = split(seqAssign, ":");
    if size(seqAssign,2) == 1, seqAssign = seqAssign'; end %if one mass 

    % detectors names are concatenated with sequence name, eg "H3S1"
    activeCollectors = extractBefore(seqAssign(:,2), seqName)';
    
    % assign massID to OP sequence table
    nMasses = size(seqAssign, 1);
    for iMass = 1:nMasses
        OPTable.(activeCollectors(iMass))(seqName) = seqAssign(iMass,1);
    end

    % record approx isotopic masses (N values) in mass table for deltas
    for iMass = 1:nMasses
        massName = extract(seqAssign(iMass,1),  lettersPattern) + ...
        extract(seqAssign(iMass,1),  digitsPattern);
        OPMasses.(activeCollectors(iMass))(seqName) = double(mass.(massName));
    end % for iMass
    
end % for iOP

OPMasses = convertvars(OPMasses, FaraNames, 'double'); % convert to double


% calculate collector deltas (amu differences between collector positions)
AxialPositionDetectorIndex = find(FaraNames == method.axialPositionDetector.Code);

OPMassMatrix = table2array(OPMasses);
nMassesPerSequence = sum(~isnan(OPMassMatrix),2);
nMassPairs = sum(nMassesPerSequence.*(nMassesPerSequence-1)/2); % total pairs
massDiffMatrix = zeros(nMassPairs, nFara-1);
massDiffColumnIdcs = [1:AxialPositionDetectorIndex-1 0 AxialPositionDetectorIndex:nFara-1];
massDiffVector = zeros(nMassPairs,1);

massDiffRow = 0;
for iSeq = 1:nSeq

    OPMassColIdcs = find(~isnan(OPMassMatrix(iSeq,:)));
    for iStartIndex = 1 : (nMassesPerSequence(iSeq)-1)
        for jEndIndex = (iStartIndex+1) : nMassesPerSequence(iSeq)
            
            OPMassStartIndex = OPMassColIdcs(iStartIndex);
            OPMassEndIndex   = OPMassColIdcs(jEndIndex);

            massDiffRow = massDiffRow + 1;
            
            if OPMassStartIndex ~= AxialPositionDetectorIndex
                massDiffMatrix(massDiffRow,massDiffColumnIdcs(OPMassStartIndex)) = -1;
            end % if iStartIndex is not axial position
            if OPMassEndIndex ~= AxialPositionDetectorIndex
                massDiffMatrix(massDiffRow,massDiffColumnIdcs(OPMassEndIndex)) = 1;
            end % if jStartIndex is not axial position

            massDiff = OPMassMatrix(iSeq,OPMassEndIndex) - ...
                       OPMassMatrix(iSeq,OPMassStartIndex);
            massDiffVector(massDiffRow) = massDiff;

        end % for end mass index            

    end % for start mass index

end % for each sequence

warning('off', 'MATLAB:rankDeficientMatrix') % suppress warning
collectorDeltasFromAxPos = (massDiffMatrix\massDiffVector)';
warning('on', 'all'); % restore warning
detectorDeltas = [collectorDeltasFromAxPos(1:AxialPositionDetectorIndex-1) ...
                   0 collectorDeltasFromAxPos(AxialPositionDetectorIndex:end)];


% create a F_ind matrix for further data reduction
% first, find unique MassIDs (stings, label for unique isotopes)
OPTableStack = stack(OPTable, FaraNames);
MassIDs = table2array(unique(OPTableStack(:,2)));
MassIDs(MassIDs == "") = []; % delete blanks
% next, create F_ind with indexes to those MassIDs
F_ind = zeros(nSeq, nFara);
for iFara = 1:nFara
    FaraColumn = OPTable.(FaraNames(iFara));

    for iMassID = 1:length(MassIDs)

        F_ind(:,iFara) = F_ind(:,iFara) + (FaraColumn == MassIDs(iMassID)) * iMassID;

    end % for iMassID

end % for iMassID

% handle baselines if present
if isfield(method, 'baselines') % if baselines present

    nBL = size(method.baselines, 2);
    BLnames = [method.baselines.Name]';
    BLTable = table('Size', [nBL, nFara], ...
        'VariableTypes', repelem("double", nFara), ...
        'VariableNames', FaraNames, ...
        'RowNames', BLnames);

    for iBL = 1:nBL

        AxMass = double(string(method.baselines(iBL).Info(4).Value)) + ...
                 double(string(method.baselines(iBL).Info(10).Value));
        BLTable{iBL,FaraNames} = AxMass + detectorDeltas;

    end % for iBL

end

% save 
method.OPTable = OPTable;
method.OPMasses = OPMasses;
method.F_ind = F_ind;
method.BLTable = BLTable;
method.MassIDs = MassIDs;
method.detectorDeltas = detectorDeltas;

end % function processMethod

