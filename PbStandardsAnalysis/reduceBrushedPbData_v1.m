function runs = reduceBrushedPbData_v1(obj, event_opj, runs, MSmethod)
% Receive handoff from raw data plotting/brushing function on UI button press.
% Record brushing results, then analyze remaining data

%% 1. Record brushing results
n.runs = size(runs, 1); % number of runs
n.BIratios = arrayfun(@(x) size(x.BIdata, 1), runs); % number of BI-ed ratios

skips = struct('PbStds', 0);
for irun1 = 1:n.runs
    
    skips(irun1).PbStds = get(runs(irun1).hPv(1), 'BrushData');
        
end
clear irun1


