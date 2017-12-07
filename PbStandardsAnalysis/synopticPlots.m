function synopticPlots(obj, event_opj) %#ok<INUSD>
session = evalin('base', 'session');
runs = evalin('base', 'runs');
uiparts = evalin('base', 'uiparts');

n.runs = size(runs,1);

close(findall(0, 'type', 'figure', 'name', 'Synoptic Plots'))

uiparts.synopticFigure = figure('Position', [100 100 1400 750], 'Name', 'Synoptic Plots');
tgroup = uitabgroup('Parent', uiparts.synopticFigure);
tab = gobjects(2, 1);

    tab(1) = uitab(tgroup, 'Title', 'Ratio vs. Intensity');
    axes('parent', tab(1))

for irun = 1:n.runs

    switch runs(irun).standard
        case 'NBS982', std = 1;
        case 'NBS981', std = 2;
    end
    
subplot(3,2,3-std); hold on
plot(runs(irun).dataDT0(runs(irun).ratioStartCycles,3), runs(irun).BIdata(:,1), '.')

subplot(3,2,5-std); hold on
plot(runs(irun).dataDT0(runs(irun).ratioStartCycles,3), runs(irun).BIdata(:,3), '.')

subplot(3,2,7-std); hold on
plot(runs(irun).dataDT0(runs(irun).ratioStartCycles,3), runs(irun).BIdata(:,4), '.')

end



