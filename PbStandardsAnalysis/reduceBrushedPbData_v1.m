function runs = reduceBrushedPbData_v1(obj, event_opj, runs, MSmethod) %#ok<INUSL>
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

stnd.logmassratios = stnd.massPb206 ./ [stnd.massPb204 stnd.massPb207 stnd.massPb208];


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

dt.in = 27.5; % nanoseconds, initial value
options = optimset('TolX', 10^-8);
[dt.best, dt.misfit] = fminsearch(@(dt) ...
                        optimize982deadtime(dt, stnd, runs, MSmethod), dt.in, options);

disp(['best fit deadtime: ' num2str(dt.best, '%2.2f') ' ns'])

% perform dead time correction with best dead time, do beam interpolation
for irun2 = 1:n.runs
    
    runs(irun2).datadt = runs(irun2).dataDT0 ./ (1 - (dt.best*1e-9)*runs(irun2).dataDT0);
    
    switch MSmethod.BImethod
        case 'Dodson'
            runs(irun2).BIdt = DodsonBI_v1(runs(irun2).datadt, MSmethod);
        case 'Quadrift'
            runs(irun2).BIdt = QuadDriftCorr_v1(runs(irun2).datadt, MSmethod);
    end % swtich beam interpolation method
    
end


%% 3.  Plot results
% alpha plots for each run (981 and 982 for now)
%figure('Name', 'Results', 'WindowSytle', 'docked');
%figure('Name', 'Results', 'WindowStyle', 'docked', 'Position', [60 60 1000 800])

hAv = gobjects(n.runs,4); % graphics array of axes handles
hPv = gobjects(n.runs,5); % graphics array of plot handles

f = figure('Position', [100 100 1400 750], 'Name', 'Results');
tgroup = uitabgroup('Parent', f);
tab = zeros(n.runs, 1);
for irun3 = 1:n.runs
    
    tab(irun3) = uitab(tgroup, 'Title', runs(irun3).name);
    axes('parent', tab(irun3))
    hold on
    
    % calculate betas - exponential fractionation factors
    
    n.ratios = size(runs(irun3).BIdt, 1);
    switch runs(irun3).standard
        case 'NBS981'
            stnd.logratios = log([stnd.nbs981r46 stnd.nbs981r76 stnd.nbs981r86]);
        case 'NBS982'
            stnd.logratios = log([stnd.nbs982r46 stnd.nbs982r76 stnd.nbs982r86]);
    end % switch case
    
    runs(irun3).beta = (repmat(stnd.logratios, n.ratios, 1) ...
                                            - log(runs(irun3).BIdt(:,[1 3 4]))) ./ ...
                        repmat(stnd.logmassratios, n.ratios, 1);
    
    % beta-46 vs. beta-86
    hAv(irun3,1) = subplot(2,3,1); hold on
    hPv(irun3,1) = plot(runs(irun3).beta(:,3), runs(irun3).beta(:,1), '.r');
    hAv(irun3,1).PlotBoxAspectRatio = [1 1 1];
    xlim = hAv(irun3,1).XLim; ylim = hAv(irun3,1).YLim;
    plot(0, 0, '+k', 'MarkerSize', 12)
    hAv(irun3,1).XLim = [min([xlim ylim]) max([xlim ylim])];
    hAv(irun3,1).YLim = hAv(irun3,1).XLim;
    line(hAv(irun3,1).XLim, hAv(irun3,1).YLim, 'LineWidth', 1, 'Color', 'k')
    
    % beta-76 vs. beta-86
    hAv(irun3,2) = subplot(2,3,4); hold on
    hPv(irun3,2) = plot(runs(irun3).beta(:,3), runs(irun3).beta(:,2), '.r');
    hAv(irun3,2).PlotBoxAspectRatio = [1 1 1];
    xlim = hAv(irun3,2).XLim; ylim = hAv(irun3,2).YLim;
    plot(0, 0, '+k', 'MarkerSize', 12)
    hAv(irun3,2).XLim = [min([xlim ylim]) max([xlim ylim])];
    hAv(irun3,2).YLim = hAv(irun3,2).XLim;
    line(hAv(irun3,2).XLim, hAv(irun3,2).YLim, 'LineWidth', 1, 'Color', 'k')
    
    % beta vs. time
    hAv(irun3,3) = subplot(2,3,[2,3]);
    hold on
    hPv(irun3, 3:5) = plot(runs(irun3).ratioTimes, runs(irun3).beta);
    legend({'\beta 204/206', '\beta 207/206', '\beta 208/206'}, 'FontSize', 12)
    
    

    
end % for irun3

tgroup.SelectedTab = tab(1);
set(findall(gcf,'type','axes'),'FontSize',16);
