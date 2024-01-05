function writeSynDataCollector(setup, filename)
%WRITESYNDATACOLLECTOR Write collector block for synthetic data file
%   Copy collector block format from Version 2.0.11 of Isolynx, 2022
%
% Written for Tripoli and CollectorRelativeEfficiencies on 5-Apr-2023
% by Noah McLean

if ~setup.writeData
    return
end

massSpec = setup.massSpec;

nIonCounters = sum(isIonCounter(massSpec.collectorArray.collectors));
nFaradays    = sum(   isFaraday(massSpec.collectorArray.collectors));
nCollectors = nIonCounters + nFaradays;
collBlock = strings(nCollectors, 6);

% ion counters
collBlock(1:nIonCounters,1)     = massSpec.ionCounterNames';
collBlock(1:nIonCounters,2)     = massSpec.ionCounterTypes';
collBlock(1:nIonCounters,3)     = compose("%1.0e", 1e11);
collBlock(1:nIonCounters,4:5)   = compose("%1.9f", 1); % gains & efficiences = 1
collBlock(1:nIonCounters,6)     = compose("%1.4f", massSpec.ionCounterDeadTimes');

% faradays
fStart = nIonCounters + 1; % start index for faradays
collBlock(fStart:end,1)   = massSpec.faradayNames';
collBlock(fStart:end,2)   = "F";
collBlock(fStart:end,3)   = compose("%1.0e", massSpec.amplifierResistance');
collBlock(fStart:end,4:5) = compose("%1.9f", 1);
collBlock(fStart:end,6)   = compose("%1.4f", 0);

collBlock = strtrim(collBlock); % remove leading and trailing whitespace
writematrix(collBlock, "../syndata/"+filename, ...
    "QuoteStrings", "none", "WriteMode", "append")

% blank line the lazy way
writematrix("", "../syndata/"+filename, ...
    "QuoteStrings", "none", "WriteMode", "append")


end

