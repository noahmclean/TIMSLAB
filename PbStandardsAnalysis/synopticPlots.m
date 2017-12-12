function synopticPlots(obj, event_opj) %#ok<INUSD>

% import data/ui structures from base workspace
session = evalin('base', 'session');
runs = evalin('base', 'runs');
uiparts = evalin('base', 'uiparts');

n.runs = size(runs,1);
cmap = colormap('lines'); % to distinguish data from different runs
lwidth = 2; % line width for uncertainty bars

% create uifigure and graphic objects structure for tabs
close(findall(0, 'type', 'figure', 'name', 'Synoptic Plots'))
uiparts.synopticFigure = figure('Position', [100 100 1400 750], 'Name', 'Synoptic Plots');
tgroup = uitabgroup('Parent', uiparts.synopticFigure);
tab = gobjects(3, 1); % 3 = ratio vs. intensity, fractionation vs. time, frac vs. temp.

% store handles for axes and plot data
hPm = gobjects(n.runs, 4); % six for first uitab
hAm = gobjects(3, 6); % likewise
hLv = gobjects(n.runs,2); % for legend

%% Ratio vs. intensity figure

tab(1) = uitab(tgroup, 'Title', 'Ratio vs. Intensity');
%axes('parent', tab(1)) % establishing this seems to be important

% set up axes, all units relative
axwidth = 0.35;   % width of plot
axheight = 0.26;  % height of plot
axsepv = 0.02;    % vertical distance between plots
axseph = 0.05;    % horizontal distance between plots
axleft = 0.07;    % distance to left side of figure
axbott = 0.09;    % distance to bottom of figure

% place axes in position, first column is 981, second column is 982
hAm(1,1) = axes('NextPlot', 'add', 'Parent', tab(1), 'Position', ... % 981 204/206
                [axleft                axbott+2*axheight+2*axsepv axwidth axheight]);
hAm(1,3) = axes('NextPlot', 'add', 'Parent', tab(1), 'Position', ... % 981 207/206
                [axleft                axbott+axheight+axsepv     axwidth axheight]);
hAm(1,5) = axes('NextPlot', 'add', 'Parent', tab(1), 'Position', ... % 981 208/206
                [axleft                axbott                     axwidth axheight]);
hAm(1,2) = axes('NextPlot', 'add', 'Parent', tab(1), 'Position', ... % 982 204/206
                [axleft+axwidth+axseph axbott+2*axheight+2*axsepv axwidth axheight]);
hAm(1,4) = axes('NextPlot', 'add', 'Parent', tab(1), 'Position', ... % 982 207/206
                [axleft+axwidth+axseph axbott+axheight+axsepv     axwidth axheight]);
hAm(1,6) = axes('NextPlot', 'add', 'Parent', tab(1), 'Position', ... % 982 208/206
                [axleft+axwidth+axseph axbott                     axwidth axheight]);

% labels, titles, and display properties
set(hAm(1, [1 2 3 4]), 'XTickLabel', [])
hAm(1,5).XLabel.String = 'Intensity (kcps)';
hAm(1,6).XLabel.String = 'Intensity (kcps)';
hAm(1,1).YLabel.String = '204/206';
hAm(1,3).YLabel.String = '207/206';
hAm(1,5).YLabel.String = '208/206';

for iax = 1:6,   hAm(1,iax).XAxis.FontSize = 12; end
for iax = 1:6,   hAm(1,iax).YAxis.FontSize = 12; end
for iax = 1:2:5,   hAm(1,iax).YLabel.FontSize = 16; end
for iax = [5 6], hAm(1,iax).XLabel.FontSize = 16; end

title1 = annotation(tab(1), 'textbox'); 
title1.Position = [axleft, axbott+3*axheight+3*axsepv, ...
                  axwidth, (0.97-(axbott+3*axheight+3*axsepv))];
title1.String = 'NBS981'; 
title1.HorizontalAlignment = 'center';
title1.VerticalAlignment = 'middle';
title1.LineStyle = 'none';

title2 = annotation(tab(1), 'textbox');
title2.Position = [axleft+axwidth+axseph, axbott+3*axheight+3*axsepv, ...
                   axwidth, (0.97-(axbott+3*axheight+3*axsepv))];
title2.String = 'NBS982'; 
title2.HorizontalAlignment = 'center';
title2.VerticalAlignment = 'middle';
title2.LineStyle = 'none';
set([title1 title2], 'FontSize', 20)

% plot deadtime-uncorrected ratios vs. intensity for each run, intensity in kcps
for irun = 1:n.runs
    
    switch runs(irun).standard
        case 'NBS982', std = 1;
        case 'NBS981', std = 2;
    end
    
hPm(irun, 1) = plot(runs(irun).dataDT0(runs(irun).ratioStartCycles,5)/1e3, ...
    runs(irun).BIdata(:,1), '.', 'Parent', hAm(1, 3-std));

hPm(irun, 2) = plot(runs(irun).dataDT0(runs(irun).ratioStartCycles,5)/1e3, ...
    runs(irun).BIdata(:,3), '.', 'Parent', hAm(1, 5-std));

hPm(irun, 3) = plot(runs(irun).dataDT0(runs(irun).ratioStartCycles,5)/1e3, ...
    runs(irun).BIdata(:,4), '.', 'Parent', hAm(1, 7-std));

end

% create large markers for legend entries in hLv, hidden behind legend
axes('Position', [0.85 0.5 0.01 0.01], 'Parent', tab(1))
hold on
for irun = 1:n.runs
    hLv(irun, 1) = plot(0, 0, '.', 'MarkerSize', 20, 'Color', cmap(irun,:));
end

% create legend, position at right of figure
pos = [(axleft + 2*axwidth + 2*axseph) axbott ...
            (1 - axleft - 2*axwidth - 3*axseph) (3*axheight+2*axsepv)]; %legend position
lgd1 = legend(hLv(:,1), {runs(:).name});
set(lgd1, 'Position', pos, 'Units', 'normalized', 'FontSize', 12)


%% Fractionation vs. temperature


tab(2) = uitab(tgroup, 'Title', 'Fractionation vs. Temp');
hAm(2,1) = axes('parent', tab(2), 'Position', ...
        [axleft axbott, 2*axwidth+axseph, (1-axbott-0.01)]);

hold on
for irun = 1:n.runs
    
    switch runs(irun).standard
        case 'NBS981', markerstr = 's';
        case 'NBS982', markerstr = 'o';
    end
    
    hPm(irun, 4) = scatter(runs(irun).temp(runs(irun).ratioStartCycles), ...
        runs(irun).beta(:,3), 100, 'Marker', markerstr, ...
        'MarkerFaceColor', cmap(irun,:), 'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', 0.7);

end

xl = xlim;
xlim([1050 xl(2)]) % 1050 is minimum realistic temperature
xlabel('Temperature (C)', 'FontSize', 20)
ylabel('\beta (208Pb/206Pb)', 'FontSize', 20)


% create large markers for legend entries in hLv, hidden behind legend
axes('Position', [0.85 0.5 0.01 0.01], 'Parent', tab(2))
hold on
for irun = 1:n.runs
    
    switch runs(irun).standard
        case 'NBS981', markerstr = 's';
        case 'NBS982', markerstr = 'o';
    end
    
    hLv(irun, 2) = scatter(0, 0, 600, 'Marker', markerstr, ...
        'MarkerFaceColor', cmap(irun,:), 'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', 0.7);
end

% make legend
warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode'); % remove tab-related warning
[lgd2, hobj2, ~, ~] = legend(hLv(:,2), {runs(:).name});
set(lgd2, 'Position', [pos(1)-0.04 pos(2) 0.15 pos(4)], 'Units', 'normalized', 'FontSize', 12)
set(findobj(hobj2,'type','patch'), 'MarkerSize',10, 'FaceAlpha', 0.7)
set(findobj(hobj2,'type','Text') , 'FontSize', 12)

%% Fractionation vs. time

    tab(3) = uitab(tgroup, 'Title', 'Fractionation vs. time');
    hAm(3,1) = axes('parent', tab(3));

for irun = 1:n.runs

    switch runs(irun).standard
        case 'NBS981', markerstr = 's';
        case 'NBS982', markerstr = 'o';
    end

    hold on

    for iratio = 1:3

        line([runs(irun).time runs(irun).time], ...
            [runs(irun).meanBeta(iratio) - 2*runs(irun).sterBeta(iratio) ...
            runs(irun).meanBeta(iratio) + 2*runs(irun).sterBeta(iratio)], ...
            'Color', cmap(iratio,:), 'LineWidth', lwidth);

        hPm(irun, iratio) = plot(runs(irun).time, runs(irun).meanBeta(iratio), ...
            'Marker', markerstr, 'MarkerFaceColor', cmap(iratio,:), ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 10);

    end % for iratio

end % for irun

datetick('x', 'dd-mmm')
ylabel('\beta', 'FontSize', 20)
xlabel('Run date', 'FontSize', 20)

