function d = assembleDataVector(data,method)
%ASSEMBLEDATAVECTOR Assemble data vector in d = g(m)
%   d is actually a struct containing the data vector and relevant
%   vectors of tags. Based off LoadMSdata_synth.m from Scott Burdick.
%   For Faraday Relative Efficiencies project, 27-Jun-2022 by Noah McLean

% determine which Faradays are used in the analysis, harmonize data/method
% matrices and indices so that columns line up with columns
[~, FaradayIndicesInData] = ...
    intersect(FaraNamesFromData, FaraNamesFromMethod, 'stable');
[FaraNamesUsed, FaradayIndicesInMethod] = ...
    intersect(FaraNamesFromMethod, FaraNamesFromData, 'stable');

F_ind = method.F_ind(:,FaradayIndicesInMethod) ;
data.BLmatrix = data.BLmatrix(:,FaradayIndicesInData);
data.OPmatrix = data.OPmatrix(:,FaradayIndicesInData);



end

