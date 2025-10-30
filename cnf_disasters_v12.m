% this code modifies v11 in following ways:
% - loop over variables 
% - quantiles

%% 0. PRELIMINARIES
%------------------------------------------------------------------
clear all; close all; clc
warning off all
format short g

%% Add path to toolbox
% addpath('functions')

% Create output folders
% mkdir('figures/cnf_disasters_v12')
mkdir('figures/cnf_disasters_v12/quantiles')

%% SECTION 1: disasters INSTRUMENT VAR
%------------------------------------------------------------------
% Load data
[disasters_xlsdata, disasters_xlstext] = xlsread('data/data_affected.xlsx','Sheet1');
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
    subplot(5,3,ii)
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


%% LOOP over disaster types
disaster_types = {'TotalaffectedFlood', 'TotalaffectedStorm', 'TotalaffectedWildfire'};
% disaster_names = {'Flood', 'Storm', 'Wildfire'};

for d = 1:length(disaster_types)
    disaster_name = disaster_types{d};
    disaster = data.(disaster_name);
    clean_name = strrep(disaster_name, 'Totalaffected', '');

    %% Plot raw disaster variable
    t = datetime(2000,1,1):calmonths(1):datetime(2019,12,1);
    figure;
    bar(t, disaster)
    title(['Number of affected: ' clean_name])
    % SaveFigure(['figures/cnf_disasters_v12/EMDAT_affected_' lower(clean_name)], 2)

    %% Exclude top 5% events
    threshold_95 = quantile(disaster, 0.95);
    disaster_top5_removed = disaster;
    disaster_top5_removed(disaster > threshold_95) = 0;

    % Plot before and after removing top 5%
    figure;
    subplot(2,1,1)
    bar(t, disaster)
    title(['All data: ' clean_name])
    subplot(2,1,2)
    bar(t, disaster_top5_removed)
    title(['Top 5% removed: ' clean_name])
    SaveFigure(['figures/cnf_disasters_v12/quantiles/EMDAT_affected_' lower(clean_name) '_noTop5'], 2)
    
    % Replace disaster
    disaster = disaster_top5_removed;

    %% Normalize disaster variable (mean of non-zero = 1)
    weights_temp = zeros(size(disaster));
    weights = zeros(size(disaster));
    minimum = min(disaster(disaster > 0));
    avg = mean(disaster(disaster > 0));

    for i = 1:length(disaster)
        weights_temp(i) = disaster(i) * minimum;
    end

    for i = 1:length(disaster)
        if disaster(i) ~= 0
            weights(i) = (nnz(disaster) / sum(weights_temp)) * weights_temp(i);
        end
    end

    disaster = weights;
    data.disaster = disaster;

    %% Create shocks and controls (same as before)
    numLags = 12;
    disaster_lags = lagmatrix(disaster, 1:numLags);
    disaster_shock = [disaster, disaster_lags];
    disaster_shock = disaster_shock(numLags+1:end, :);

    numLagsCV = 12;
    lagged_vars = {'INDPRO','CPIAUCSL','DGS1'};

    for i = 1:length(lagged_vars)
        var = lagged_vars{i};
        data.([var '_lags']) = lagmatrix(data.(var), 1:numLagsCV);
        data.([var '_lags']) = data.([var '_lags'])(numLags+1:end, :);
    end

    controls = [data.INDPRO_lags, data.CPIAUCSL_lags, data.DGS1_lags];

    %% Responses
    response_vars = {'disaster','INDPRO','CPIAUCSL','DGS1'};
    for i = 1:length(response_vars)
        var = response_vars{i};
        data.(['Y' num2str(i)]) = data.(var)(numLags+1:end, :);
    end

    T = size(data.Y2, 1);
    data.C = [ones(T,1)];
    data.X1 = [data.C, disaster_shock, controls];

    %% Run LP regressions
    horizons = 46;
    dependent_indices = 1:4;
    confidence_levels = [1, 1.645];

    for idx = dependent_indices
        Y = data.(['Y' num2str(idx)]);
        X = data.X1;

        for h = 0:horizons-1
            Xh = X(2:end-h, :);
            Y_t_1 = Y(1:end-h-1);
            Y_t_h = Y(h+2:end);
            Y_diff = Y_t_h - Y_t_1;

            r = rank(Xh);
            cond_num = rcond(Xh'*Xh);

            if r < size(Xh,2) || cond_num < 1e-10
                warning(['X matrix at horizon ' num2str(h) ' is rank deficient. Using pinv.']);
                coeff = pinv(Xh) * Y_diff;
            else
                coeff = regress(Y_diff, Xh);
            end

            uhat = Y_diff - Xh * coeff;
            se = NeweyWest(uhat, Xh);

            coeff = coeff(:)';
            se = se(:)';

            data.(['coeff' num2str(idx)])(h+1,:) = coeff;
            for c = 1:length(confidence_levels)
                z = confidence_levels(c);
                data.(['coeff' num2str(idx) '_high' num2str(c)])(h+1,:) = coeff + z * se;
                data.(['coeff' num2str(idx) '_low'  num2str(c)])(h+1,:) = coeff - z * se;
            end
        end
    end

    %% Plot
    figure;
    avg_rounded_k = round(avg / 1000);  % Round to nearest thousand
    sgtitle(sprintf('Impulse: %s \nAvg. of non-zero event: %dK people', clean_name, avg_rounded_k));
    main_line_color = [0, 0, 0];
    shade_color = [0.8, 0.8, 0.8];

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
        high2 = data.([base '_high2'])(:,2)';
        low2  = data.([base '_low2'])(:,2)';
        high1 = data.([base '_high1'])(:,2)';
        low1  = data.([base '_low1'])(:,2)';
        coeff = data.(base)(:,2);

        fill([0:horizons-1, fliplr(0:horizons-1)], [high2, fliplr(low2)], ...
            shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.3); 
        fill([0:horizons-1, fliplr(0:horizons-1)], [high1, fliplr(low1)], ...
            shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        plot(0:horizons-1, coeff, 'Color', main_line_color, 'LineWidth', 2); 
        line([0, horizons-1], [0, 0], 'Color', 'black', 'LineWidth', 1);
        title(plot_configs{i}.title, 'FontSize', 16);
        set(gca, 'FontSize', 14);
        grid on;
    end

    % Save figure with disaster name
    SaveFigure(sprintf('figures/cnf_disasters_v12/quantiles/IRFs_affected_%s_noTop5', lower(clean_name)), 2);
    % SaveFigure(sprintf('figures/cnf_disasters_v12/IRFs_affected_%s', lower(clean_name)), 2);
end
