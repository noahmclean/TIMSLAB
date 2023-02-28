%% How long should I measure a baseline?

% constants
kB = 1.38064852*10^-23;  % Boltzman's constant
T = 290; % Kelvin
R = 1E11;  % resistor size
amplifierNoiseVariance = 4*kB*T*R;  %Johnson-Nyquist noise

cpsPerVolt =  6.241509074460763e+18 / R;  % elementary charge, Ohm's Law

ttotal = 100; % seconds
tints = 1000; % total time integrations
tvec = linspace(0, ttotal, tints+2)';
tvec = tvec(2:end-1); % trim first and last time (zero baseline/on-peak time)

intv = [0.025 0.050 0.100 0.200 0.400 0.800]* (R/10^11); % intensities, Volts
% coded inputs to intv are 10^11 - equivalent intensities
nints = length(intv);

% do calculations for all possible beam intensities 
unctot = zeros(tints, nints);
for ii = 1:nints
    
    varBL = amplifierNoiseVariance./tvec; % variance from baseline meas.
    varOP = (intv(ii)/cpsPerVolt + amplifierNoiseVariance)./(ttotal-tvec); % variance from on peak meas.
    unctot(:,ii) = sqrt(varBL + varOP)/intv(ii); % total relative uncertainty (1-sigma)
    
end

[mins, mindx] = min(unctot, [], 1); % find optimum (minimum uncertainty) BL vs. OP timing

%% Figure and color map
f1 = figure('Position', [10 10 1000*0.9 450*0.95], 'Units', 'pixels');

map = brighten(colormap(winter(nints+2)), 0); % choose and adjust color map
map = map(2:end-1,:); % reject more saturated colors at start and end of map


%% Plot 1 - logarithmic y-axis

ax1 = subplot(1,2,1);
hold on
p1 = zeros(nints,1); p2 = p1;
for ii = 1:nints
    p1(ii) = plot(tvec, unctot(:,ii), 'Color', map(ii,:), 'LineWidth', 2);
    p2(ii) = plot(tvec(mindx(ii)), mins(ii), '.', 'MarkerSize', 18, 'Color', map(ii,:));
end

ax1.YScale = 'log';
ax1.FontSize = 14;
ax1.XLabel.String = ['seconds of baseline'];% per ' num2str(ttotal, '%0.0f') ' seconds of total time'];
ax1.YLabel.String = 'fractional uncertainty';
%ax1.Title.String = {'Fractional Uncertainty in Baseline-Corrected Intensity for',
%    [num2str(ttotal, '%0.0f') ' Seconds of Total Analytical Time']};
ax1.LineWidth = 1;
ax1.Box = 'on';
pos = ax1.Position;
ax1.Position = pos + [0 -0.00 0 -0.04];
%ax.YLim = [10^2 10^6];

lgd1 = legend(p1, cellstr(num2str(intv', '%0.3f V')));
lgd1.FontSize = 16;
lgd1.Location = 'north';
lgd1.Box = 'on';
%lgd.Interpreter = 'latex';
%lgdpos = lgd.Position;
%lgd.Position = lgdpos .* [0.7 0.9 2 1.7];

%% Plot 2 - same info, linear y-axis

ax2 = subplot(1,2,2);
hold on
p3 = zeros(nints,1); p4 = p3;
for ii = 1:nints
    p3(ii) = plot(tvec, unctot(:,ii), 'Color', map(ii,:), 'LineWidth', 2);
    p4(ii) = plot(tvec(mindx(ii)), mins(ii), '.', 'MarkerSize', 18, 'Color', map(ii,:));
end

%ax = gca;
%ax2.YScale = 'log';
ax2.FontSize = 14;
ax2.XLabel.String = ['seconds of baseline'];% per ' num2str(ttotal, '%0.0f') ' seconds of total time'];
ax2.YLabel.String = 'fractional uncertainty';
%ax2.Title.String = {'Fractional Uncertainty in Baseline-Corrected Intensity for',
%    [num2str(ttotal, '%0.0f') ' Seconds of Total Analytical Time']};
ax2.LineWidth = 1;
ax2.Box = 'on';
ax2.YLim = [0 0.0008];
pos = ax2.Position;
ax2.Position = pos + [0 -0.00 0 -0.04];

lgd2 = legend(p3, cellstr(num2str(intv', '%0.3f V')));
lgd2.FontSize = 16;
lgd2.Location = 'north';
lgd2.Box = 'on';
%lgd.Interpreter = 'latex';
%lgdpos = lgd.Position;
%lgd.Position = lgdpos .* [0.7 0.9 2 1.7];

annotation('textbox', [0 0.9 1 0.1], ...
    'String', ['Fractional Uncertainty in Baseline-Corrected Intensity for ', ...
      num2str(ttotal, '%0.0f')  ' Seconds of Total Analytical Time, ', ...
      'R = 10^{',num2str(log10(R), '%0.0f') ,'} \Omega'], ...
    'EdgeColor', 'none', 'FontSize', 18, ...
    'HorizontalAlignment', 'center')