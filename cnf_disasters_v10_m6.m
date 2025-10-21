% this code modifies v10 in following ways:
% - remove 6 first months of observations

%% 0. PRELIMINARIES
%------------------------------------------------------------------
clear all; close all; clc
warning off all
format short g

%% Add path to toolbox

% Create output folders
mkdir('figures/cnf_disasters_v10_m6')

%% SECTION 1: disasters INSTRUMENT VAR
%------------------------------------------------------------------
% Load data
[disasters_xlsdata, disasters_xlstext] = xlsread('data/data_disasters_v10_m6.xlsx','Sheet1');
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
FigSize(26,18)
for ii=1:disasters_nvar
    subplot(4,3,ii)
    H(ii) = plot(disasters_DATA.(disasters_vnames{ii}),'LineWidth',3,'Color',cmap(1));
    title(disasters_vnames(ii)); 
    DatesPlot(disasters_datesnum(1),disasters_nobs,6,'m') % Set the x-axis label 
    grid on; 
end
clf('reset')

% Set up LP (Eickmeier, et al., 2024) 
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

% First difference 
fd_vars = [{'UNRATE','DGS1','EBP'}, log_vars]; % Include level vars + log-transformed ones

for i = 1:length(fd_vars)
    var = fd_vars{i};
    data.([var '_FD']) = diff(data.(var));
end

% Disaster variable
disaster = data.DISASTERS;

time = datetime(2000,7,1):calmonths(1):datetime(2019,12,1);

figure;
bar(time, disaster)
SaveFigure('figures/cnf_disasters_v10_m6/EMDAT',2)

% Disaster shock: contemporaneous + 9 lags 
numLags=9;
disaster_lags = lagmatrix(disaster, 1:numLags);
disaster_shock = [disaster, disaster_lags];
disaster_shock = disaster_shock(numLags+1:end, :);

% Lagged control variables
numLagsCV=3;
% List of control variables (FD)
lagged_vars = {'UNRATE_FD','CPIAUCSL_FD','HousePr_FD','StockPr_FD','VIX_FD',...
                'DGS1_FD','EBP_FD','INDPRO_FD','GFC','CCPI_FD','CPIUFDSL_FD'};

% Create lagged versions
for i = 1:length(lagged_vars)
    var = lagged_vars{i};
    data.([var '_lags']) = lagmatrix(data.(var), 1:numLagsCV);
end

% Add linear trend
T = size(disaster_shock, 1);
data.TREND = (1:T)';

% Adjust the lenght of the series 
for i = 1:length(lagged_vars)
    var = lagged_vars{i};
    lagged_var = [var '_lags'];
    
    % GFC starts one period later
    if strcmp(var, 'GFC')
        data.(lagged_var) = data.(lagged_var)(numLags+1:end, :);
    else
        data.(lagged_var) = data.(lagged_var)(numLags:end, :);
    end
end

controls = [data.UNRATE_FD_lags, data.CPIAUCSL_FD_lags, data.HousePr_FD_lags, data.StockPr_FD_lags, data.VIX_FD_lags, data.DGS1_FD_lags, data.GFC_lags];

% Variables for response (level variables)
response_vars = {'UNRATE','CPIAUCSL','HousePr','StockPr','VIX','DGS1','GFC','INDPRO','EBP','CCPI','CPIUFDSL'};

% Create response variables
for i = 1:length(response_vars)
    var = response_vars{i};
    data.(['Y' num2str(i+1)]) = data.(var)(numLags+1:end, :);
end

% Deterministic component that comprises a constant and a linear trend
T = size(data.Y2, 1);
data.C = [ones(T,1), data.TREND(end-T+1:end)];

% Create X that contains constant, control variable and shock
% Base X matrix
data.X1 = [data.C, disaster_shock, controls];

% Additional X matrices with specific lags 
data.X9 = [data.C disaster_shock controls data.INDPRO_FD_lags];
data.X10 = [data.C disaster_shock controls data.EBP_FD_lags];
data.X11 = [data.C disaster_shock controls data.CCPI_FD_lags];
data.X12 = [data.C disaster_shock controls data.CPIUFDSL_FD_lags];

horizons=36; 

%% Long-run differences
dependent_indices = 2:12;  % Corresponds to Y2 to Y12
confidence_levels = [1, 1.645];  % 68% and 90% confidence intervals

for idx = dependent_indices
    Y = data.(['Y' num2str(idx)]);

    % Select appropriate X matrix
    if ismember(idx, [9, 10, 11, 12])
        X = data.(['X' num2str(idx)]);
    else
        X = data.X1;
    end

    for h = 1:horizons
        Xh = X(1:end-h, :);
        Y_t = Y(1:end-h, :);
        Y_t_h = Y(h+1:end, :);
        Y_diff = Y_t_h - Y_t;

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

        % Ensure row vector format
        coeff = coeff(:)';
        se = se(:)';

        % Store coefficients
        data.(['coeff' num2str(idx)])(h,:) = coeff;

        for c = 1:length(confidence_levels)
            z = confidence_levels(c);
            data.(['coeff' num2str(idx) '_high' num2str(c)])(h,:) = coeff + z * se;
            data.(['coeff' num2str(idx) '_low' num2str(c)])(h,:)  = coeff - z * se;
        end
    end
end

% Plot
figure;
main_line_color = [0, 0, 0]; % black
shade_color = [0.8, 0.8, 0.8]; % light gray
% Define plot configurations
plot_configs = {
    struct('title', 'INDPRO', 'base', 'coeff9'), ...    
    struct('title', 'CPI', 'base', 'coeff3'), ...
    struct('title', 'UNRATE', 'base', 'coeff2'), ...  
    struct('title', '1Y Treasury Bond', 'base', 'coeff7')
};
for i = 1:length(plot_configs)
    subplot(2, 2, i);
    hold on;
    base = plot_configs{i}.base;
    title_text = plot_configs{i}.title;
    % Dynamically access fields
    high2 = data.([base '_high2'])(:,3)';
    low2  = data.([base '_low2'])(:,3)';
    high1 = data.([base '_high1'])(:,3)';
    low1  = data.([base '_low1'])(:,3)';
    coeff = data.(base)(:,3);
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
SaveFigure('figures/cnf_disasters_v10_m2/IRFs_disasters', 2);
clf('reset')