function method = processMethod(method, FaraNames, collectorDeltas)
%PROCESSMETHOD Extract useful information from parsed method structure
%
%   Then append the useful info to the method structure

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

OPMasses = convertvars(OPMasses, FaraNames, 'double'); % conver to double

% % NOTE: need to calculate mass differences in a different way for 
% non-consecutive masses.  Maybe make a distance/mileage chart between all
% collector pairs?
% % calculate mass differences
% OPMasses = double(table2array(OPMasses));
% OPMassDifTable = zeros(size(OPMasses) - [0 1]);
% for iSeq = 1:nSeq
%     seqMass = OPMasses(iSeq,:);
%     massPresent = seqMass > 0;
%     massDifPlace = massPresent(2:end) & massPresent(1:end-1);
%     massVec = seqMass(massPresent);
%     massDifs = massVec(2:end)-massVec(1:end-1);
%     OPMassDifTable(iSeq, massDifPlace) = massDifs;
% end % for iMass
% 
% OPMassDif = zeros(1, nFara-1);
% for iFara = 1:nFara-1
% 
%     massDifFara = OPMassDifTable(:,iFara);
%     massDifFara = massDifFara(massDifFara > 0);
%     if ~all(massDifFara == massDifFara(1))
%         disp('problem with OPMassDifTable')
%     else 
%         OPMassDif(iFara) = massDifFara(1);
%     end
% 
% end % for iFara

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
        BLTable{iBL,FaraNames} = AxMass + collectorDeltas;

    end % for iBL

end


method.OPTable = OPTable;
method.OPMasses = OPMasses;
method.F_ind = F_ind;
method.BLTable = BLTable;
method.MassIDs = MassIDs;

end % function processMethod

