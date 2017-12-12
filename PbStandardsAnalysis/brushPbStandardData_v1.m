function runs = brushPbStandardData_v1(runs, MSmethod)
%  Interactively discard bad cycles of Pb data
%  Record these in a field in the runs structure for each run, called skipv
%  skipv = 0 for good data to use in calculations, 1 for bad data to be excluded.

n.runs = size(runs, 1);


%% Set up some data brushing

clear gcf

hAv = gobjects(n.runs,8); % graphics array of axes handles
hPv = gobjects(n.runs,8); % graphics array of plot handles


close(findall(0, 'Type', 'figure', 'Name', 'Raw Data'))
uiparts.rawDataFigure = figure('Position', [100 100 1400 750], 'Name', 'Raw Data');

tgroup = uitabgroup('Parent', uiparts.rawDataFigure, 'TabLocation', 'top');
tab = zeros(n.runs, 1);

for irun = 1:n.runs
    
    tab(irun) = uitab(tgroup, 'Title', runs(irun).name, 'Tag', num2str(irun));
    axes('parent', tab(irun))
    
    % Determine when BI ratios start by enumerating the first cycle of each ratio
    
    % 1. Count blocks and cycles
    numBlocks = size(runs(irun).dataRaw,1)/MSmethod.cyclesPerBlock; %total (plus fractional) blocks of data
    nBlocks = floor(numBlocks); % complete blocks of data
    partialBlockCycles = max(rem(size(runs(irun).dataRaw,1),MSmethod.cyclesPerBlock)-1,0);
    % 2. make a vector of the first cycle of each BI pair for full blocks
    cycleCount = repmat(1:(MSmethod.cyclesPerBlock-1),1,nBlocks);
    blockTotal = repelem(0:MSmethod.cyclesPerBlock:(MSmethod.cyclesPerBlock*nBlocks-1), MSmethod.cyclesPerBlock - 1);
    ratioStartCycles = cycleCount + blockTotal; % cycle indices for first cycle in BI ratios
    % 3. append cycles from partial blocks to end; does nothing for partialBlockCycles = 0.
    ratioStartCycles = [ratioStartCycles  (1:partialBlockCycles) + MSmethod.cyclesPerBlock*nBlocks]; %#ok<AGROW>

    runs(irun).ratioTimes = runs(irun).secs(ratioStartCycles);
    runs(irun).ratioStartCycles = ratioStartCycles; % record this for later, where it's useful.
    
    hAv(irun,1) = subplot(4,2,1); hPv(irun,1) = plot(runs(irun).ratioTimes/60, runs(irun).BIdata(:,1), '.'); ylabel('204/206')
    hAv(irun,2) = subplot(4,2,3); hPv(irun,2) = plot(runs(irun).ratioTimes/60, runs(irun).BIdata(:,2), '.'); ylabel('205/206')
    hAv(irun,3) = subplot(4,2,5); hPv(irun,3) = plot(runs(irun).ratioTimes/60, runs(irun).BIdata(:,3), '.'); ylabel('207/206')
    hAv(irun,4) = subplot(4,2,7); hPv(irun,4) = plot(runs(irun).ratioTimes/60, runs(irun).BIdata(:,4), '.'); ylabel('208/206')

    hAv(irun,5) = subplot(4,2,2); hPv(irun,5) = plot(runs(irun).ratioTimes/60, runs(irun).temp(ratioStartCycles), '.')     ; ylabel('temp')
    hAv(irun,6) = subplot(4,2,4); hPv(irun,6) = plot(runs(irun).ratioTimes/60, runs(irun).dataRaw(ratioStartCycles,1), '.'); ylabel('204 cps')
    hAv(irun,7) = subplot(4,2,6); hPv(irun,7) = plot(runs(irun).ratioTimes/60, runs(irun).dataRaw(ratioStartCycles,2), '.'); ylabel('205 cps')
    hAv(irun,8) = subplot(4,2,8); hPv(irun,8) = plot(runs(irun).ratioTimes/60, runs(irun).dataRaw(ratioStartCycles,5), '.'); ylabel('208 cps')
    
    set(hPv(irun,:), 'Tag', num2str(irun));
    
    linkaxes([hAv(irun,1) hAv(irun,2) hAv(irun,3) hAv(irun,4) hAv(irun,5) hAv(irun,6) hAv(irun,7) hAv(irun,8)], 'x');
    
    set(hAv(irun, 1:8), 'XLim', [min(runs(irun).ratioTimes/60) max(runs(irun).ratioTimes/60)])
    
    temp.xdatastr = ['runs(' int2str(irun) ').ratioTimes'];
    set(hPv(irun,1:8), 'XDataSource', temp.xdatastr)
    
    temp.h = zoom;
    set(temp.h, 'Motion', 'vertical', 'Enable', 'on');
    brush('on')
    set(brush, 'ActionPostCallback', {@uponBrushing, hPv})
    
end % plotting for irun = 1:nruns

% Transfer hPv to the 'runs' structure
for irun = 1:n.runs
    runs(irun).hPv = hPv(irun,:);
end

close(findall(0, 'Type', 'figure', 'Name', 'Controls'))
uiparts.ControlButton = figure('name', 'Controls', ...
                                'WindowStyle', 'normal', 'Position', [50 50 200 200]);
uicontrol('Style', 'pushbutton', 'String', 'Reduce Data', ...
                                 'Position', [10 60 100 20], ... 
                                 'Callback', {@reduceBrushedPbData_v1, runs, MSmethod});
uicontrol('Style', 'pushbutton', 'String', 'Synoptic Plots', ...
                                 'Position', [10 10 100 20], ...
                                 'Callback', {@synopticPlots});

                             
% keep uiparts in base workspace for other UI functions
uiparts.rawDataTgroup = tgroup;
assignin('base', 'uiparts', uiparts);


%%%%%%%%%%%%%%%%%%%%%%%%
%%% Callbacks %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

function uponBrushing(obj, event_opj, hPv)
% obj: handle to the figure the has been clicked
% event_obj: object conting struct of event data

justBrushed = findobj(event_opj.Axes, 'Type', 'line');
currentBrushing = get(justBrushed, 'BrushData');
brushedRun = justBrushed.Parent.Parent.Tag; % from 'Tag' for tab
allLinesInRun = findobj(hPv, 'Tag', brushedRun);
set(allLinesInRun, 'BrushData', currentBrushing);


end % function: link plot data brushing callback


end % function
