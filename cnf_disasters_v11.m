% this code modifies v10 in following ways:
% - variables in levels
% - limit to 3 basic variables and see the effect of each disaster
% - horizon from 0 (!)

%% 0. PRELIMINARIES
%------------------------------------------------------------------
clear all; close all; clc
warning off all
format short g

%% Add path to toolbox
% addpath('functions')

% Create output folders
mkdir('figures/cnf_disasters_v11')

%% SECTION 1: disasters INSTRUMENT VAR
%------------------------------------------------------------------
% Load data
[disasters_xlsdata, disasters_xlstext] = xlsread('data/data_affected_storm.xlsx','Sheet1');
disasters_dates = disasters_xlstext(3:end,1);
disasters_datesnum = Date2Num(disasters_dates, 'm');
disasters_vnames_long = disasters_xlstext(1,2:end);
disasters_vnames = disasters_xlstext(2,2:end);
disasters_nvar = length(disasters_vnames);
disasters_data = Num2NaN(disasters_xlsdata);

% Build structure
for ii = 1:disasters_nvar
    disasters_DATA.(disasters_vnames{ii}) = disasters_data(:,ii);
end
disasters_nobs = size(disasters_data,1);

% Plot
figure;
for ii=1:disasters_nvar
    subplot(4,3,ii)
    H(ii) = plot(disasters_DATA.(disasters_vnames{ii}),'LineWidth',3,'Color',cmap(1));
    title(disasters_vnames(ii)); 
    DatesPlot(disasters_datesnum(1),disasters_nobs,6,'m') % Set the x-axis label 
    grid on; 
end

% Initialize structure
data = struct();

% Assign each variable to the structure (Rename)
for i = 1:disasters_nvar
    disasters_varname = disasters_vnames{i};
    data.(disasters_varname) = disasters_xlsdata(:,i);
end

% Log level
log_vars = {'CPIAUCSL','HousePr','StockPr','VIX','INDPRO','CCPI', 'CPIUFDSL'};

for i = 1:length(log_vars)
    disasters_var = log_vars{i};
    data.(disasters_var) = 100 * log(data.(disasters_var));
end

% Disaster variable
% disaster = data.Flood;
% disaster = data.Storm;
% disaster = data.Wildfire;
% disaster = data.TotalaffectedFlood;
disaster = data.TotalaffectedStorm;
% disaster = data.TotalaffectedWildfire;

t = datetime(2000,1,1):calmonths(1):datetime(2019,12,1);

figure;
bar(t, disaster)
SaveFigure('figures/cnf_disasters_v11/EMDAT_affected_storm',2)

%% Exclude outliers (Interquartile Range):
% P5 = quantile(disaster, 0.05);
% P95 = quantile(disaster, 0.95);
% IQR = P95 - P5;
% 
% lower_bound = P5 - 1.5 * IQR;
% upper_bound = P95 + 1.5 * IQR;
% 
% % Logical mask for inliers
% inlier_idx = (disaster >= lower_bound) & (disaster <= upper_bound);
% 
% % Copy the original data
% disasterr = disaster;
% 
% % Replace outliers with 0
% disasterr(~inlier_idx) = 0;
% 
% figure;
% subplot(2,1,1)
% bar(t, disaster)
% title('all data')
% 
% subplot(2,1,2)
% bar(t, disasterr)
% title('Without Outliers')
% SaveFigure('figures/cnf_disasters_v11/EMDAT_affected_storm',2)
% 
% % Run the regression without outliers
% disaster = disasterr;

%% Normalize disaster variable so the average of non-zero shocks equals 1 (event/affected)
weights_temp = zeros(size(disaster));
weights = zeros(size(disaster));
minimum = min(disaster(disaster > 0));
avg = mean(disaster(disaster>0));

for i = 1:length(disaster)
    weights_temp(i) = disaster(i) * minimum;
end

for i = 1:length(disaster)
    if disaster(i) ~= 0
        weights(i) = (nnz(disaster) / sum(weights_temp)) * weights_temp(i);
    end
end
disaster = weights;
data.disaster=disaster;

% Disaster shock: contemporaneous + 9 lags 
numLags=12;
disaster_lags = lagmatrix(disaster, 1:numLags);
disaster_shock = [disaster, disaster_lags];
disaster_shock = disaster_shock(numLags+1:end, :);

% Lagged control variables
numLagsCV=12;
% List of control variables (FD)
lagged_vars = {'INDPRO','CPIAUCSL','DGS1'};

% Create lagged versions
for i = 1:length(lagged_vars)
    var = lagged_vars{i};
    data.([var '_lags']) = lagmatrix(data.(var), 1:numLagsCV);
end

% Adjust the lenght of the series 
for i = 1:length(lagged_vars)
    var = lagged_vars{i};
    lagged_var = [var '_lags'];
    data.(lagged_var) = data.(lagged_var)(numLags+1:end, :);
end

controls = [data.INDPRO_lags, data.CPIAUCSL_lags, data.DGS1_lags];

% Variables for response (level variables)
response_vars = {'disaster','INDPRO','CPIAUCSL','DGS1'};

% Create response variables
for i = 1:length(response_vars)
    var = response_vars{i};
    data.(['Y' num2str(i)]) = data.(var)(numLags+1:end, :);
end

% Deterministic component that comprises a constant and a linear trend
T = size(data.Y2, 1);
data.C = [ones(T,1)];

% Create X that contains constant, control variable and shock
% Base X matrix
data.X1 = [data.C, disaster_shock, controls];

horizons=46; 

%% Long-run differences
dependent_indices = 1:4;  % Corresponds to Y1 to Y4
confidence_levels = [1, 1.645];  % 68% and 90% confidence intervals

for idx = dependent_indices
    Y = data.(['Y' num2str(idx)]);

    % Select appropriate X matrix
    X = data.X1;
    for h = 0:horizons-1
        Xh = X(2:end-h, :);         % X_t aligned with Y_t_1
        Y_t_1 = Y(1:end-h-1);       % Y_{t-1}
        Y_t_h = Y(h+2:end);         % Y_{t+h}
        Y_diff = Y_t_h - Y_t_1;

        % Check rank and condition number
        r = rank(Xh);
        cond_num = rcond(Xh'*Xh);

        if r < size(Xh,2) || cond_num < 1e-10
            warning(['Matrix X' num2str(idx) ' at horizon ' num2str(h) ' is rank deficient or ill-conditioned. Using pinv.']);
            coeff = pinv(Xh) * Y_diff;
        else
            coeff = regress(Y_diff, Xh);
        end

        uhat = Y_diff - Xh * coeff;
        se = NeweyWest(uhat, Xh);

        % Ensure row vectors
        coeff = coeff(:)';
        se = se(:)';

        % Store coefficients
        data.(['coeff' num2str(idx)])(h+1,:) = coeff;
        for c = 1:length(confidence_levels)
            z = confidence_levels(c);
            data.(['coeff' num2str(idx) '_high' num2str(c)])(h+1,:) = coeff + z * se;
            data.(['coeff' num2str(idx) '_low'  num2str(c)])(h+1,:) = coeff - z * se;
        end
    end
end

% Plot 1
figure;
sgtitle(sprintf('Impulse: \nAvg. affected: %.0f people', avg))
main_line_color = [0, 0, 0]; % black
shade_color = [0.8, 0.8, 0.8]; % light gray
% Define plot configurations
plot_configs = {
    struct('title', 'Affected', 'base', 'coeff1'), ...
    struct('title', 'INDPRO', 'base', 'coeff2'), ... 
    struct('title', 'CPI', 'base', 'coeff3'), ...
    struct('title', '1Y Treasury Bond', 'base', 'coeff4')
};
for i = 1:length(plot_configs)
    subplot(2, 2, i);
    hold on;
    base = plot_configs{i}.base;
    title_text = plot_configs{i}.title;
    % Dynamically access fields
    high2 = data.([base '_high2'])(:,2)';
    low2  = data.([base '_low2'])(:,2)';
    high1 = data.([base '_high1'])(:,2)';
    low1  = data.([base '_low1'])(:,2)';
    coeff = data.(base)(:,2);
    % Plot shaded areas and main line
    fill([0:horizons-1, fliplr(0:horizons-1)], [high2, fliplr(low2)], ...
        shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.3); 
    fill([0:horizons-1, fliplr(0:horizons-1)], [high1, fliplr(low1)], ...
        shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    plot(0:horizons-1, coeff, 'Color', main_line_color, 'LineWidth', 2); 
    line([0, horizons-1], [0, 0], 'Color', 'black', 'LineWidth', 1);
    title(title_text, 'FontSize', 20);
    set(gca, 'FontSize', 16);
    grid on;
end
SaveFigure('figures/cnf_disasters_v11/IRFs_affected_storm', 2);