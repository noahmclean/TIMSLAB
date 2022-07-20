function [d, data] = assembleDataVector(data,method)
%ASSEMBLEDATAVECTOR Assemble data vector in d = g(m)
%   d is actually a struct containing the data vector and relevant
%   vectors of tags. Based off LoadMSdata_synth.m from Scott Burdick.
%   For Faraday Relative Efficiencies project, 27-Jun-2022 by Noah McLean

% determine which Faradays are used in the analysis, harmonize data/method
% matrices and indices so that columns line up with columns
DetNamesFromMethod = string(method.OPTable.Properties.VariableNames);
DetNamesFromData = data.collectorNames;
[~, FaradayIndicesInData] = ...
    intersect(DetNamesFromData, DetNamesFromMethod, 'stable');
[FaraNamesUsed, FaradayIndicesInMethod] = ...
    intersect(DetNamesFromMethod, DetNamesFromData, 'stable');

BLTable = method.BLTable{:,FaradayIndicesInMethod};
F_ind = method.F_ind(:,FaradayIndicesInMethod) ;
data.BLmatrix = data.BLmatrix(:,FaradayIndicesInData);
data.OPmatrix = data.OPmatrix(:,FaradayIndicesInData);

nOPseq = size(F_ind,1);
nBLseq = size(BLTable,1);
nDet = size(F_ind,2);

% determine baseline references
BLdetIsRef = zeros(size(BLTable)); % true if baseline int is referenced
for iOPseq = 1:nOPseq
    
    % which baselines are referenced by each sequence?
    refString = string(method.onpeaks(iOPseq).Info(12).Value);    
    %BLrefs is a vector with the BL squence indices referenced in iOPseq
    BLrefs = double(extractAfter(split(refString, ", "), "BL"));
    
    % make a flag for the BLTable with which intensities are used
    detInSeq = F_ind(iOPseq,:) > 0;
    BLdetIsRef(BLrefs, :) = BLdetIsRef(BLrefs, :) | detInSeq;

end
OPdetIsRef = F_ind > 0; % true of OP int is referenced


%% make masks for BL and OP data used in inversion

BLdataIsRef = logical(BLdetIsRef(data.BLSeqIdx,:));
OPdataIsRef = logical(OPdetIsRef(data.OPSeqIdx,:));



%% correct bad 2-second relay settle time for gains on select samples
% goodGains contains 6-second RST ccgains measured on 08-Jul-2022 (mean of 100)

if any([ "SmEfficiency_Bead3Run2-1809.TXT" "SmEfficiency_Bead3Run3-1810.TXT"] ...
               == string(data.header.fileName) )
goodGains = [0.9964579, 1.0118078, 0.9885882, 1.0086657, 1.0000000, ...
             1.0149499, 0.9886156, 0.9841766, 0.9649764];


% undo bad gains, apply good gains
data.BLmatrix = data.BLmatrix ./ (data.FaradayGains');
data.BLmatrix = data.BLmatrix .* goodGains;
data.OPmatrix = data.OPmatrix ./ (data.FaradayGains');
data.OPmatrix = data.OPmatrix .* goodGains;

end % if any Sm were measured with bad 2-second RST ccgains


%% assemble baseline data into data vector/struct

% initialize d -- structure with data (intensities)
% and metadata needed to evaluate model to calculate dhat
d.int   = []; % intensity
d.time  = []; % time
d.det   = []; % detector index
d.mass  = []; % isotope/species mass 
d.iso   = []; % isotope index (for OP)
d.block = []; % block index
d.isOP  = []; % is data point on OP measurement?

for iDet = 1:nDet 

    detRefs = BLdataIsRef(:,iDet); % where is this detector refd in BL?
    nrefs   = sum(detRefs);
    d.int   = [d.int; data.BLmatrix(detRefs, iDet)];
    d.time  = [d.time; data.BLtime(detRefs)];
    d.det   = [d.det; iDet*ones(nrefs,1)];
    BLseqs  = data.BLSeqIdx(detRefs); % BL seqs referenced
    d.mass  = [d.mass; BLTable( BLseqs, iDet )];
    d.iso   = [d.iso; zeros(nrefs,1)];
    d.block = [d.block; data.BLserial(detRefs, 1)];
    d.isOP  = [d.isOP; false(nrefs,1)];

end


%% append on-peak data onto d vector/struct

for iDet = 1:nDet 

    detRefs = OPdataIsRef(:,iDet); % where is this detector refd in F_ind?
    nrefs   = sum(detRefs);
    d.int   = [d.int; data.OPmatrix(detRefs, iDet)];
    d.time  = [d.time; data.OPtime(detRefs)];
    d.det   = [d.det; iDet*ones(nrefs,1)];
    OPseqs  = data.OPSeqIdx(detRefs); % OP seqs referenced
    d.mass  = [d.mass; method.OPMasses{ OPseqs, iDet }];
    d.iso   = [d.iso; F_ind(OPseqs,iDet)];
    d.block = [d.block; data.OPserial(detRefs, 1)];
    d.isOP  = [d.isOP; true(nrefs,1)];

end



end % function