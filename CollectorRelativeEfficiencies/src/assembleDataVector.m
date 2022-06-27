function d = assembleDataVector(data,method)
%ASSEMBLEDATAVECTOR Assemble data vector in d = g(m)
%   d is actually a struct containing the data vector and relevant
%   vectors of tags. Based off LoadMSdata_synth.m from Scott Burdick.
%   For Faraday Relative Efficiencies project, 27-Jun-2022 by Noah McLean

% determine which Faradays are used in the analysis
FaraNamesFromMethod = string(method.OPTable.Properties.VariableNames);
FaraNamesFromData = data.collectorNames;
[FaraNamesUsed, FaradayIndicesInData] = intersect(FaraNamesFromMethod, FaraNamesFromData);

FarsUsed = find(sum(method.F_ind));    
F_ind = method.F_ind(:,FarsUsed) ;




end

