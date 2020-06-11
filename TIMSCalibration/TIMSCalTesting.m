%% Call and test functions for TIMS Calibration
% Master file for new functions, written starting 6/2020
% Analyze gains and standard data, including parsing functions
% for eventual data processing (Tripoli, Redux) functionality

%% Section 1: Gains

% Data folder of .raw files 
dataFolderPath = ...
'/Volumes/groups/GEOL/Mass_Spec/new_backup/';
%'/Users/noahmc/Google Drive/SPOCK-536 Tracer Calibration Data/';

% find CC Gains files in fol
userpath(dataFolderPath);
gainsFolderStruct = dir([dataFolderPath 'CCGains*.raw']);

nGains = size(gainsFolderStruct,1);

% query each .raw folder for its gains
deleteGainsList = [];
gains = struct('data', {}, 'mean', {}, 'lstd', {}, 'lste', {}); %preallocate
for iFolder = 1:nGains
    
    % determine gains folder, read data
    gainsFolderPath = [dataFolderPath gainsFolderStruct(iFolder).name '/'];
    
    gainsFileInfo = dir([gainsFolderPath '*.xls']);
    
    if ~size(gainsFileInfo,1) % if there's no Excel file in the folder
        nGains = nGains - 1;
        deleteGainsList = [deleteGainsList iFolder]; %#ok<AGROW>
        continue % skip data compilation
    end
    
    disp(num2str(iFolder))
    
    gainsRawData = readmatrix([gainsFolderPath gainsFileInfo.name], ...
                           'Sheet', 'Gains', 'Range', 'E20:N119');
    
    % clean up empty files and NaNs
    if all(all(isnan(gainsRawData))) % if no data -- all NaNs
        nGains = nGains - 1;
        deleteGainsList = [deleteGainsList iFolder]; %#ok<AGROW>
        continue
    end
    
    gainsData = gainsRawData(:,[1:4 6:10]); % delete IC column
    gainsData(any(isnan(gainsData), 2), :) = []; % delete rows w/ NaNs
    
    gains(iFolder).data = gainsData;
    gains(iFolder).mean = geomean(gainsData);
    gains(iFolder).lstd = std(log(gainsData));
    gains(iFolder).lste = gains(iFolder).lstd/sqrt(size(gainsData,1));
    gains(iFolder).date = gainsFileInfo.datenum;

end % 

% clean up gains runs identified as incomplete
gainsFolderStruct(deleteGainsList) = [];
gains(deleteGainsList) = [];


%% Pick out and sort data for plotting

% sort gains structure by dateData
dateData = vertcat(gains.date);
[dateData, idx] = sort(dateData);
gains = gains(idx);
gainsFolderStruct = gainsFolderStruct(idx);

% delete gain 16 for now...
gains(16) = [];
dateData(16) = [];

% pick out gains statistics
meanData = vertcat(gains.mean);
lstdData = vertcat(gains.lstd);
lsteData = vertcat(gains.lste);

% express data as ppm for easy comparison on plots
meanDataPPM = (meanData - mean(meanData))./mean(meanData) * 1e6;
lsteDataPPM = lsteData * 1e6; % standard error as ppm, linearized approx


%% Make some quick plots

% gains to plot: High and Low collectors (relative to Axial = 1)
g2plot = [1:4 6:9];
gainLabels = {'L5/Ax', 'L4/Ax', 'L3/Ax', 'L2/Ax', 'H1/Ax', 'H2/Ax', 'H3/Ax', 'H4/Ax'};

% x-axis is 1:n or by date (can cluster)
%xData = 1:nGains;
xData = dateData';


% set up a subplot-type figure but with adjacent axes
topbuf = 0.02;
botbuf = 0.05;
totvbuf = topbuf + botbuf;
plotHeight = (1-totvbuf)/8;

% Vertical stack of gains
set(0, 'units', 'pixels')
monitorSize = get(0, 'screensize');
figure('Position', [0 0 monitorSize(3)/2 monitorSize(4)], ...
       'Units', 'pixels', 'DefaultAxesFontSize', 16);

for iPlot = 1:8
    
    ax(iPlot) = axes('Position', [0.1 botbuf+(iPlot-1)*plotHeight 0.85 plotHeight], ...
                     'Units', 'normalized');
    ylim([-100 100]);
    set(gca, 'XTickLabel', []) % remove xtick labels (restore bottom plot's below)
    ylabel(gainLabels(iPlot), 'FontSize', 18)
    
    hold on
    plot(xData, meanDataPPM(:,g2plot(iPlot)), '-', 'LineWidth', 2) % gains
    line([xData; xData], ... % 2SE uncertainties
         [meanDataPPM(:,g2plot(iPlot))' - 2*lsteDataPPM(:,g2plot(iPlot))'; 
          meanDataPPM(:,g2plot(iPlot))' + 2*lsteDataPPM(:,g2plot(iPlot))'], ...
          'LineWidth', 2, 'Color', 'k');
    
end

% restore bottom axis ticks
set(ax(1), 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
datetick(ax(1))

% reset xlimits for plots because relabelling can rescale them
for iPlot = 2:8
    ax(iPlot).XLim = ax(1).XLim;
    ax(iPlot).XTick = ax(1).XTick;
end

% remove duplicate axis labels
for iPlot = 1:7
    ax(iPlot).YTick = ax(iPlot).YTick(1:end-1);
end


