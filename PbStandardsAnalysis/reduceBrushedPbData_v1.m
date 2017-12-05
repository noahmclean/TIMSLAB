function runs = reduceBrushedPbData_v1(obj, event_opj, runs, MSmethod) %#ok<INUSL>
% Receive handoff from raw data plotting/brushing function on UI button press.
% Record brushing results, then analyze remaining data

% standard and physical constant values:
% from Condon/McLean spike calibration papers
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

stnd.logmassratios = log(stnd.massPb206 ./ [stnd.massPb204 stnd.massPb207 stnd.massPb208]);


%% 1. Record brushing results
n.runs = size(runs, 1); % number of runs
n.BIratios = arrayfun(@(x) size(x.BIdata, 1), runs); % number of BI-ed ratios

for irun = 1:n.runs
    
    runs(irun).skips = logical(get(runs(irun).hPv(1), 'BrushData'))';
        
end % for irun, record brushing results


%% 2. Do deadtime calculations
% First pass: do all calculations with NBS982, then correct NBS 981 data as 
% a secondary standard

dt.init = 27.5; % nanoseconds, initial value
options = optimset('TolX', 10^-8);
[dt.best, dt.misfit] = fminsearch(@(dt) ...
                        optimize982deadtime(dt, stnd, runs, MSmethod), dt.init, options);

disp(['best fit deadtime: ' num2str(dt.best, '%2.2f') ' ns'])

% perform dead time correction with best dead time, do beam interpolation
for irun = 1:n.runs
    
    runs(irun).datadt = runs(irun).dataDT0 ./ (1 - (dt.best*1e-9)*runs(irun).dataDT0);
    
    switch MSmethod.BImethod
        case 'Dodson'
            runs(irun).BIdt = DodsonBI_v1(runs(irun).datadt, MSmethod);
        case 'Quadrift'
            runs(irun).BIdt = QuadDriftCorr_v1(runs(irun).datadt, MSmethod);
    end % swtich beam interpolation method
    
end % for irun, do dead time calculations


%% Some fractionation calculations

    % calculate alphas - linear fractionation factors, and
    %           betas - exponential fractionation factors
for irun = 1:n.runs
    
    n.ratios = size(runs(irun).BIdt, 1);
    switch runs(irun).standard
        case 'NBS981'
            stnd.ratios    = [stnd.nbs981r46 stnd.nbs981r76 stnd.nbs981r86];
            stnd.logratios = log([stnd.nbs981r46 stnd.nbs981r76 stnd.nbs981r86]);
        case 'NBS982'
            stnd.ratios    = [stnd.nbs982r46 stnd.nbs982r76 stnd.nbs982r86];
            stnd.logratios = log([stnd.nbs982r46 stnd.nbs982r76 stnd.nbs982r86]);
    end % switch case
    
    runs(irun).alpha = repmat([-1/2 1 1/2], n.ratios, 1) .* ...
                      (repmat(stnd.ratios,  n.ratios, 1) ./ ...
                       runs(irun).BIdt(:, [1 3 4]) - 1 );
    
    runs(irun).beta = (repmat(stnd.logratios, n.ratios, 1) ...
        - log(runs(irun).BIdt(:,[1 3 4]))) ./ ...
        repmat(stnd.logmassratios, n.ratios, 1);
    
    if ~isempty(runs(irun).skips) % if the user brushed some data
        runs(irun).meanAlpha = mean(runs(irun).alpha(~runs(irun).skips,:));
        runs(irun).sterAlpha  = std(runs(irun).alpha(~runs(irun).skips,:))/...
                                sqrt(length(runs(irun).alpha(~runs(irun).skips,:)));
        runs(irun).meanBeta = mean(runs(irun).beta(~runs(irun).skips,:));
        runs(irun).sterBeta  = std(runs(irun).beta(~runs(irun).skips,:))/...
                               sqrt(length(runs(irun).beta(~runs(irun).skips,:)));
        ratioTemps = runs(irun).temp(runs(irun).ratioStartCycles); % just the temps at ratio times
        runs(irun).meanTemp = mean(ratioTemps(~runs(irun).skips)); % skips is tied to ratios not cycles
    else
        runs(irun).meanAlpha = mean(runs(irun).alpha);
        runs(irun).sterAlpha  = std(runs(irun).alpha)/...
                                sqrt(length(runs(irun).alpha));
        runs(irun).meanBeta = mean(runs(irun).beta);
        runs(irun).sterBeta  = std(runs(irun).beta)/...
                               sqrt(length(runs(irun).beta));
        runs(irun).meanTemp = mean(runs(irun).temp);
    end % if user brushed some data
        
end % for irun, fractionation calculations

% calculate session weighed mean fractionation values
session.wtdMeanAlpha   = sum(cat(1, runs.meanAlpha) ./ cat(1,runs.sterAlpha).^2) ./ ...
                         sum(                     1 ./ cat(1,runs.sterAlpha).^2);
session.wtdMeanAlpha1s = sqrt(1./sum(1./cat(1, runs.sterAlpha).^2));
session.redChiSqAlpha  = sum((cat(1,runs.meanAlpha)-repmat(session.wtdMeanAlpha1s,n.runs,1)).^2 ./ ...
                              cat(1,runs.sterAlpha).^2 ) / (n.runs - 1);

session.wtdMeanBeta =   sum(cat(1, runs.meanBeta)   ./ cat(1,runs.sterBeta).^2) ./ ...
                        sum(                      1 ./ cat(1,runs.sterBeta).^2);
session.wtdMeanBeta1s  = sqrt(1./sum(1./cat(1, runs.sterBeta).^2));
session.redChiSqBeta  = sum((cat(1,runs.meanBeta)-repmat(session.wtdMeanBeta1s,n.runs,1)).^2 ./ ...
                              cat(1,runs.sterBeta).^2 ) / (n.runs - 1);
session.meanTemp = mean(cat(1, runs.meanTemp));

%% 3.  Plot results
% alpha plots for each run (981 and 982 for now)
%figure('Name', 'Results', 'WindowSytle', 'docked');
%figure('Name', 'Results', 'WindowStyle', 'docked', 'Position', [60 60 1000 800])

hAv = gobjects(n.runs,4); % graphics array of axes handles
hPv = gobjects(n.runs,5); % graphics array of plot handles

% close previous 
close(findall(0, 'type', 'figure', 'name', 'Run Results'))

uiparts.runResultsFigure = figure('Position', [100 100 1400 750], 'Name', 'Run Results');
tgroup = uitabgroup('Parent', uiparts.runResultsFigure);
tab = gobjects(n.runs, 1);

%setup for data table on each tab:
alphastr = '%1.2f'; % format string for alphas
betastr  = '%1.2f'; % format string for betas

str6  = num2str(session.wtdMeanAlpha(1)*100, alphastr);
str7  = num2str(session.wtdMeanAlpha1s(1)*200, alphastr);
str8  = num2str(session.wtdMeanBeta(1), betastr);
str9  = num2str(session.wtdMeanBeta1s(1)*2, betastr);
str10 = num2str(session.redChiSqBeta(1), '%4.1f');
str16  = num2str(session.wtdMeanAlpha(2)*100, alphastr);
str17  = num2str(session.wtdMeanAlpha1s(2)*200, alphastr);
str18  = num2str(session.wtdMeanBeta(2), betastr);
str19  = num2str(session.wtdMeanBeta1s(2)*2, betastr);
str20  = num2str(session.redChiSqBeta(2), '%4.1f');
str26  = num2str(session.wtdMeanAlpha(3)*100, alphastr);
str27  = num2str(session.wtdMeanAlpha1s(3)*200, alphastr);
str28  = num2str(session.wtdMeanBeta(3), betastr);
str29  = num2str(session.wtdMeanBeta1s(3)*2, betastr);
str30  = num2str(session.redChiSqBeta(3), '%4.1f');
strT2  = num2str(session.meanTemp, '%4.0f');

for irun = 1:n.runs
    
    tab(irun) = uitab(tgroup, 'Title', runs(irun).name);
    axes('parent', tab(irun))
    hold on
    
    % beta-46 vs. beta-86
    hAv(irun,1) = subplot(2,3,1); hold on
    hPv(irun,1) = plot(runs(irun).beta(:,3), runs(irun).beta(:,1), '.r');
    hAv(irun,1).PlotBoxAspectRatio = [1 1 1];
    xlim = hAv(irun,1).XLim; ylim = hAv(irun,1).YLim;
    plot(0, 0, '+k', 'MarkerSize', 12)
    hAv(irun,1).XLim = [min([xlim ylim]) max([xlim ylim])];
    hAv(irun,1).YLim = hAv(irun,1).XLim;
    line(hAv(irun,1).XLim, hAv(irun,1).YLim, 'LineWidth', 1, 'Color', 'k')
    xlabel('\beta 208/206'); ylabel('\beta 204/206');
    
    % beta-76 vs. beta-86
    hAv(irun,2) = subplot(2,3,4); hold on
    hPv(irun,2) = plot(runs(irun).beta(:,3), runs(irun).beta(:,2), '.r');
    hAv(irun,2).PlotBoxAspectRatio = [1 1 1];
    xlim = hAv(irun,2).XLim; ylim = hAv(irun,2).YLim;
    plot(0, 0, '+k', 'MarkerSize', 12)
    hAv(irun,2).XLim = [min([xlim ylim]) max([xlim ylim])];
    hAv(irun,2).YLim = hAv(irun,2).XLim;
    line(hAv(irun,2).XLim, hAv(irun,2).YLim, 'LineWidth', 1, 'Color', 'k')
    xlabel('\beta 208/206'); ylabel('\beta 207/206');
    
    % beta vs. time
    hAv(irun,3) = subplot(2,3,[2,3]);
    hold on
    hPv(irun, 3:5) = plot(runs(irun).ratioTimes/60, runs(irun).beta);
    legend({'\beta 204/206', '\beta 207/206', '\beta 208/206'}, 'FontSize', 12)
    xlabel('run time (minutes)')
    ylabel('\beta')
    

    
    str1  = num2str(runs(irun).meanAlpha(1)*100, alphastr);
    str2  = num2str(runs(irun).sterAlpha(1)*200, alphastr);
    str3  = num2str(runs(irun).meanBeta(1),      betastr);
    str4  = num2str(runs(irun).sterBeta(1)*2,    betastr);
    str5  = num2str(0);
    str11  = num2str(runs(irun).meanAlpha(2)*100, alphastr);
    str12  = num2str(runs(irun).sterAlpha(2)*200, alphastr);
    str13  = num2str(runs(irun).meanBeta(2),      betastr);
    str14  = num2str(runs(irun).sterBeta(2)*2,    betastr);
    str15  = num2str(0);
    str21  = num2str(runs(irun).meanAlpha(3)*100, alphastr);
    str22  = num2str(runs(irun).sterAlpha(3)*200, alphastr);
    str23  = num2str(runs(irun).meanBeta(3),      betastr);
    str24  = num2str(runs(irun).sterBeta(3)*2,    betastr);
    str25  = num2str(0);
    strT1  = num2str(runs(irun).meanTemp, '%4.0f');
    
    % build a long string for annotation
    tbl = '\begin{tabular}{l r@{$\pm$}l r@{$\pm$}l l c r@{$\pm$}l r@{$\pm$}l l } ';
    tbl = [tbl '  & \multicolumn{5}{c}{this run} & & \multicolumn{5}{c}{grand mean} \\ '];
    tbl = [tbl '& \multicolumn{2}{c}{$\alpha$} & \multicolumn{2}{c}{$\beta$} & $\chi^2_{red}$ '];
    tbl = [tbl '& & \multicolumn{2}{c}{$\alpha$} & \multicolumn{2}{c}{$\beta$} & $\chi^2_{red}$ \\'];
    tbl = [tbl ' \cline{2-6}\cline{8-12} '];
    tbl = [tbl '204/206   &'  str1 '&'  str2 '&'  str3 '&'  str4 '&'  str5 '& &' ];
    tbl = [tbl      str6 '&'  str7 '&'  str8 '&'  str9 '&' str10 '\\ '];
    tbl = [tbl '207/206   &' str11 '&' str12 '&' str13 '&' str14 '&' str15 '& &' ];
    tbl = [tbl     str16 '&' str17 '&' str18 '&' str19 '&' str20 '\\ '];
    tbl = [tbl '208/206   &' str21 '&' str22 '&' str23 '&' str24 '&' str25 '& &' ];
    tbl = [tbl     str26 '&' str27 '&' str28 '&' str29 '&' str30 '\\ \\ '];
    tbl = [tbl 'Temp. & \multicolumn{2}{l}{' strT1 '} & \multicolumn{2}{c}{ } & & & \multicolumn{2}{l}{' strT2 '} '];
    tbl = [tbl '\end{tabular}'];
    annotation(tab(irun), 'textbox', [0.34 0.152 0.55 0.3], 'String', tbl, ...
         'FontSize', 18, 'BackgroundColor', [0.999 0.999 0.999], 'Interpreter', 'latex')


end % for irun, plot results

tgroup.SelectedTab = tab(1);
set(findall(gcf,'type','axes'),'FontSize',16);

assignin('base', 'runs', runs);
assignin('base', 'session', session);
assignin('base', 'uiparts', uiparts);
