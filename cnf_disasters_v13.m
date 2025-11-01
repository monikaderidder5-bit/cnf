% this code modifies v12 in following ways:
% - independent variables in FD (change in horizons)
% - dependent variable is long differences Y_{t+h} - Y_{t-1}

%% 0. PRELIMINARIES
%------------------------------------------------------------------
clear all; close all; clc
warning off all
format short g

%% Add path to toolbox
% addpath('functions')

% Create output folders
mkdir('figures/cnf_disasters_v13')

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

%% LOOP over disaster types
disaster_types = {'TotalaffectedFlood', 'TotalaffectedStorm', 'TotalaffectedWildfire'};

for d = 1:length(disaster_types)
    disaster_name = disaster_types{d};
    disaster = data.(disaster_name);
    clean_name = strrep(disaster_name, 'Totalaffected', '');

    %% Plot raw disaster variable
    t = datetime(2000,1,1):calmonths(1):datetime(2019,12,1);
    figure;
    bar(t, disaster)
    title(['Number of affected: ' clean_name])
    SaveFigure(['figures/cnf_disasters_v13/EMDAT_affected_' lower(clean_name)], 2)

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
    disaster_shock = disaster_shock(numLags+2:end, :);

    numLagsCV = 12;
    lagged_vars = {'INDPRO_FD','CPIAUCSL_FD','DGS1_FD'};

    for i = 1:length(lagged_vars)
        var = lagged_vars{i};
        data.([var '_lags']) = lagmatrix(data.(var), 1:numLagsCV);
        data.([var '_lags']) = data.([var '_lags'])(numLags+1:end, :);
    end

    controls = [data.INDPRO_FD_lags, data.CPIAUCSL_FD_lags, data.DGS1_FD_lags];

    %% Responses
    response_vars = {'disaster','INDPRO','CPIAUCSL','DGS1'};
    for i = 1:length(response_vars)
        var = response_vars{i};
        data.(['Y' num2str(i)]) = data.(var)(numLags+1:end-1, :);
    end

    T = size(controls, 1);
    data.C = [ones(T,1)];
    data.X1 = [data.C, disaster_shock, controls];

    %% Run LP regressions
    horizons = 46;
    dependent_indices = 1:4;
    confidence_levels = [1, 1.645];

    for idx = dependent_indices
        Y = data.(['Y' num2str(idx)]);
        X = data.X1;

        for h = 1:horizons
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

            data.(['coeff' num2str(idx)])(h,:) = coeff;
            for c = 1:length(confidence_levels)
                z = confidence_levels(c);
                data.(['coeff' num2str(idx) '_high' num2str(c)])(h,:) = coeff + z * se;
                data.(['coeff' num2str(idx) '_low'  num2str(c)])(h,:) = coeff - z * se;
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
    % saveas(d,sprintf(['figures/cnf_disasters_v13' num2str(d) '/affected.png'],d));
    SaveFigure(sprintf('figures/cnf_disasters_v13/IRFs_affected_%s', lower(clean_name)), 2);
end
