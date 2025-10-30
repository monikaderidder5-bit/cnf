% this code modifies v12 in following ways:
% - deaths

%% 0. PRELIMINARIES
%------------------------------------------------------------------
clear all; close all; clc
warning off all
format short g

%% Add path to toolbox
% addpath('functions')

% Create output folders
% mkdir('figures/cnf_disasters_v12')
mkdir('figures/cnf_disasters_v12/quantiles/deaths')

%% SECTION 1: disasters INSTRUMENT VAR
%------------------------------------------------------------------
% Load data
[disasters_xlsdata, disasters_xlstext] = xlsread('data/data_deaths.xlsx','Sheet1');
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
disaster_types = {'TotaldeathsFlood', 'TotaldeathsStorm', 'TotaldeathsWildfire'};
% disaster_names = {'Flood', 'Storm', 'Wildfire'};

for d = 1:length(disaster_types)
    disaster_name = disaster_types{d};
    disaster = data.(disaster_name);
    clean_name = strrep(disaster_name, 'Totaldeaths', '');

    % %% Plot raw disaster variable
    t = datetime(2000,1,1):calmonths(1):datetime(2019,12,1);
    % figure;
    % bar(t, disaster)
    % title(['Number of deaths: ' clean_name])
    % SaveFigure(['figures/cnf_disasters_v12/EMDAT_deaths_' lower(clean_name)], 2)

    %% Define quantiles to keep (top X% most deadly events)
    quantile_levels = [0.99, 0.95, 0.80, 0.50, 0.01];  % Corresponding to 1%, 5%, 20%, 50%, 99%

    for q = 1:length(quantile_levels)
        q_level = quantile_levels(q);
        q_str = strrep(num2str(round((1 - q_level) * 100)), '.', '_');  % e.g., "1", "5"
    
        %% === FILTER: Keep only top X% of disasters ===
        threshold = quantile(disaster, q_level);
        top_idx = (disaster >= threshold);
    
        disaster_q = zeros(size(disaster));
        disaster_q(top_idx) = disaster(top_idx);
    
        %% Plot raw vs top-X% data
        figure;
        subplot(2,1,1)
        bar(t, disaster)
        title(['All data: ' clean_name])
    
        subplot(2,1,2)
        bar(t, disaster_q)
        title([sprintf('Top %.0f%% most deadly: ', (1 - q_level)*100), clean_name])
    
        SaveFigure(sprintf('figures/cnf_disasters_v12/quantiles/deaths/EMDAT_deaths_%s_top%s', ...
            lower(clean_name), q_str), 2)
    
        %% Normalize filtered disaster
        weights_temp = zeros(size(disaster_q));
        weights = zeros(size(disaster_q));
        if any(disaster_q > 0)
            minimum = min(disaster_q(disaster_q > 0));
            avg = mean(disaster_q(disaster_q > 0));
            for i = 1:length(disaster_q)
                weights_temp(i) = disaster_q(i) * minimum;
            end
            for i = 1:length(disaster_q)
                if disaster_q(i) ~= 0
                    weights(i) = (nnz(disaster_q) / sum(weights_temp)) * weights_temp(i);
                end
            end
        end
        disaster_q = weights;
        data.disaster = disaster_q;

    %% Create shocks and controls (same as before)
    numLags = 12;
    disaster_lags = lagmatrix(disaster_q, 1:numLags);
    disaster_shock = [disaster_q, disaster_lags];
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
    avg_rounded_k = round(avg);  % Avg. deaths (in K)
    sgtitle(sprintf('Impulse: %s \nTop %.0f%% events only | Avg. deaths: %d people', ...
        clean_name, (1 - q_level)*100, avg_rounded_k));
    main_line_color = [0, 0, 0];
    shade_color = [0.8, 0.8, 0.8];

    plot_configs = {
        struct('title', 'deaths', 'base', 'coeff1'), ...
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

    SaveFigure(sprintf('figures/cnf_disasters_v12/quantiles/deaths/IRFs_deaths_%s_top%s', ...
        lower(clean_name), q_str), 2);
    end
end
