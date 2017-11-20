function runs = reduceBrushedPbData_v1(obj, event_opj, runs, MSmethod)
% Receive handoff from raw data plotting/brushing function on UI button press.
% Record brushing results, then analyze remaining data

% standard and physical constant values:
stnd.nbs981r46 = 0.0590074;
stnd.nbs981r76 = 0.914683;
stnd.nbs981r86 = 2.1681;

stnd.nbs982r46 = 0.0272058;
stnd.nbs982r76 = 0.466967;
stnd.nbs982r86 = 1.000249;

stnd.massPb204 = 203.9730436;
stnd.massPb205 = 204.9744818;
stnd.massPb206 = 205.9744653;
stnd.massPb207 = 206.9758969;
stnd.massPb208 = 207.9766521;

%% 1. Record brushing results
n.runs = size(runs, 1); % number of runs
n.BIratios = arrayfun(@(x) size(x.BIdata, 1), runs); % number of BI-ed ratios

for irun1 = 1:n.runs
    
    runs(irun1).skips = logical(get(runs(irun1).hPv(1), 'BrushData'))';
        
end
clear irun1


%% 2. Do deadtime calculations
% First pass: do all calculations with NBS982, then correct NBS 981 data as 
% a secondary standard

dt.in = 25; % nanoseconds, initial value
options = optimset('TolX', 10^-8);
[dt.best, dt.misfit] = fminsearch(@(dt) ...
                        optimize982deadtime(dt, stnd, runs, MSmethod), dt.in, options);

disp(['best fit deadtime: ' num2str(dt.best, '%2.2f') ' ns'])

