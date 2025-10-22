% this code modifies v9 in following ways:
% - 1 unit impact of storms, floods and extreme tempretures
% - Based on eye-balling of Hack et al (2023) that defence spending shock 
% can mute out the FFR responce. Maybe, it is like that in the reality,
% that strong fiscal responce in that case do not cause inflationary effects. 

%% 0. PRELIMINARIES
%------------------------------------------------------------------
clear all; close all; clc
warning off all
format short g

%% Add path to toolbox

% Create output folders
mkdir('figures/cnf_disasters_v10')

%% SECTION 1: disasters INSTRUMENT VAR
%------------------------------------------------------------------
% Load data
[disasters_xlsdata, disasters_xlstext] = xlsread('data/data_disasters_v10.xlsx','Sheet1');
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

time = datetime(2000,1,1):calmonths(1):datetime(2019,12,1);

figure;
bar(time, disaster)
SaveFigure('figures/cnf_disasters_v10/EMDAT',2)

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
dependent_indices = 2:12;  % Corresponds to Y2 to Y11
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

% Plot 1
figure;
main_line_color = [0, 0, 0]; % black
shade_color = [0.8, 0.8, 0.8]; % light gray
% Define plot configurations
plot_configs = {
    struct('title', 'UNRATE', 'base', 'coeff2'), ...    
    struct('title', 'CPI', 'base', 'coeff3'), ...
    struct('title', 'House Price', 'base', 'coeff4'), ...
    struct('title', 'Stock Price', 'base', 'coeff5'), ...
    struct('title', 'VIX', 'base', 'coeff6'), ...    
    struct('title', '1Y Treasury Bond', 'base', 'coeff7')
};
for i = 1:length(plot_configs)
    subplot(2, 3, i);
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
SaveFigure('figures/cnf_disasters_v10/IRFs_disasters1', 2);
clf('reset')

% Plot 2
figure;
main_line_color = [0, 0, 0]; % black
shade_color = [0.8, 0.8, 0.8]; % light gray
% Define plot configurations
plot_configs = {
    struct('title', 'INDPRO', 'base', 'coeff9'), ...
    struct('title', 'EBP', 'base', 'coeff10'), ...
    struct('title', 'Commodity Price', 'base', 'coeff11'), ...
    struct('title', 'CPI Food', 'base', 'coeff12')
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
SaveFigure('figures/cnf_disasters_v10/IRFs_disasters2', 2);
clf('reset')

%% SECTION 2: MONETARY POLICY (MP) INSTRUMENT FOR Short Rate
%------------------------------------------------------------------
[MP_xlsdata, MP_xlstext] = xlsread('data/data_JAR_u1.xlsx','Sheet1');
MP_dates = MP_xlstext(3:end,1);
MP_datesnum = Date2Num(MP_dates, 'm');
MP_vnames_long = MP_xlstext(1,2:end);
MP_vnames = MP_xlstext(2,2:end);
MP_nvar = length(MP_vnames);
MP_data = Num2NaN(MP_xlsdata);

% Structure
for ii = 1:MP_nvar
    MP_DATA.(MP_vnames{ii}) = MP_data(:,ii);
end
MP_nobs = size(MP_data,1);

% Log level
log_vars = {'CPIAUCSL','HousePr','StockPr','VIX','INDPRO','CCPI', 'CPIUFDSL'};

for i = 1:length(log_vars)
    MP_var = log_vars{i};
    MP_DATA.(MP_var) = 100 * log(MP_DATA.(MP_var));
end

% VAR setup
MP_Xvnames_long = {'1Y Treasury Bond';'INDPRO';'CPI';'UNRATE';'House Price';'Stock Price';'VIX';'GFC';'EBP';'Commodity Price';'CPI Food'};
MP_Xvnames = {'DGS1';'INDPRO';'CPIAUCSL';'UNRATE';'HousePr';'StockPr';'VIX';'GFC';'EBP';'CCPI';'CPIUFDSL'};
MP_Xnvar = length(MP_Xvnames);
MP_X = nan(MP_nobs,MP_Xnvar);
for ii=1:MP_Xnvar
    MP_X(:,ii) = MP_DATA.(MP_Xvnames{ii});
end
MP_IVvnames = {'u1'};
MP_IV = MP_DATA.(MP_IVvnames{1});

% Estimate VAR
MP_det = 1;
MP_nlags = 12;
[MP_VAR, MP_VARopt] = VARmodel(MP_X,MP_nlags,MP_det);
MP_VARopt.vnames = MP_Xvnames_long;
MP_VARopt.nsteps = 36;
MP_VARopt.impact    = 1;
MP_VARopt.quality = 2;
MP_VARopt.FigSize = [26,12];
MP_VARopt.firstdate = MP_datesnum(1);
MP_VARopt.frequency = 'm';
MP_VAR.IV = MP_IV;
MP_VARopt.ident = 'iv';
MP_VARopt.snames = {'MP shock','','','','','','','','','','',''};
MP_VARopt.method = 'wild';
MP_VARopt.pctg   = 90;

% Compute IRFs
[MP_IR, MP_VAR] = VARir(MP_VAR, MP_VARopt);
[MP_IRinf, MP_IRsup, MP_IRmed, MP_IRbar] = VARirband(MP_VAR, MP_VARopt);
VARirplot(MP_IRbar, MP_VARopt, MP_IRinf, MP_IRsup);

% Plot IRFs
cmap = lines;
FigSize(26, 12);
SwatheOpt = PlotSwatheOption;
SwatheOpt.swathecol = [0.85 0.85 0.85];
SwatheOpt.linecol = [0.2 0.2 0.2];
for ii = 1:MP_Xnvar
    subplot(4, 3, ii);
    PlotSwathe(MP_IRbar(:,ii), [MP_IRinf(:,ii) MP_IRsup(:,ii)], SwatheOpt); 
    hold on;
    plot(1:MP_VARopt.nsteps, zeros(MP_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
    plot(1:MP_VARopt.nsteps, MP_IRbar(:,ii), 'Color', cmap(1,:), 'LineWidth', 2);
    title([MP_Xvnames_long{ii} ' to ' MP_VARopt.snames{1}], 'FontWeight', 'bold', 'FontSize', 10); 
    xlim([1 MP_VARopt.nsteps]);
    grid on;
    set(gca, 'Layer', 'top');
end
SaveFigure('figures/cnf_disasters_v10/IRFs_u1', 2);
clf('reset')

%% Counterfactuals
% Regress CPI Food of disasters_IRbar on CPI Food of MP_IRbar
beta = - (MP_IRbar(:,11) \ data.coeff12(:,3));  % negative to minimize combination
% You get: combined = disasters + beta * MP ≈ 0
reg = data.coeff12(:,3) + beta * MP_IRbar(:,11);
plot(reg)

% Plot IRFs
FigSize(26, 12);
% Assume cmap, SwatheOpt, MP_IRbar, MP_IRinf, MP_IRsup, MP_VARopt are already defined
% CPI Food
PlotSwathe(data.coeff12(:,3), [data.coeff12_low1(:,3)  data.coeff12_high1(:,3)], SwatheOpt); 
hold on;
plot(1:MP_VARopt.nsteps, zeros(MP_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
plot(1:MP_VARopt.nsteps, data.coeff12(:,3), 'Color', cmap(1,:), 'LineWidth', 2);
plot(1:MP_VARopt.nsteps, reg , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');  % red dashed difference
title([MP_Xvnames_long{11} ' to disaster shock'], 'FontWeight', 'bold', 'FontSize', 10); 
xlim([1 MP_VARopt.nsteps]);
SaveFigure('figures/cnf_disasters_v10/cnf_u1', 2);
clf('reset')

%% SECTION 3: MONETARY POLICY (MP) INSTRUMENT VAR FOR FG
%------------------------------------------------------------------
[MP2_xlsdata, MP2_xlstext] = xlsread('data/data_JAR_u2.xlsx','Sheet1');
MP2_dates = MP2_xlstext(3:end,1);
MP2_datesnum = Date2Num(MP2_dates, 'm');
MP2_vnames_long = MP2_xlstext(1,2:end);
MP2_vnames = MP2_xlstext(2,2:end);
MP2_nvar = length(MP2_vnames);
MP2_data = Num2NaN(MP2_xlsdata);

% Structure
for ii = 1:MP2_nvar
    MP2_DATA.(MP2_vnames{ii}) = MP2_data(:,ii);
end
MP2_nobs = size(MP2_data,1);

% Log level
log_vars = {'CPIAUCSL','HousePr','StockPr','VIX','INDPRO','CCPI', 'CPIUFDSL'};

for i = 1:length(log_vars)
    MP2_var = log_vars{i};
    MP2_DATA.(MP2_var) = 100 * log(MP2_DATA.(MP2_var));
end

% VAR setup
MP2_Xvnames_long = {'1Y Treasury Bond';'INDPRO';'CPI';'UNRATE';'House Price';'Stock Price';'VIX';'GFC';'EBP';'Commodity Price';'CPI Food'};
MP2_Xvnames = {'DGS1';'INDPRO';'CPIAUCSL';'UNRATE';'HousePr';'StockPr';'VIX';'GFC';'EBP';'CCPI';'CPIUFDSL'};
MP2_Xnvar = length(MP2_Xvnames);
MP2_X = nan(MP2_nobs,MP2_Xnvar);
for ii=1:MP2_Xnvar
    MP2_X(:,ii) = MP2_DATA.(MP2_Xvnames{ii});
end
MP2_IVvnames = {'u2'};
MP2_IV = MP2_DATA.(MP2_IVvnames{1});

% Estimate VAR
MP2_det = 1;
MP2_nlags = 12;
[MP2_VAR, MP2_VARopt] = VARmodel(MP2_X,MP2_nlags,MP2_det);
MP2_VARopt.vnames = MP2_Xvnames_long;
MP2_VARopt.nsteps = 36;
MP2_VARopt.impact    = 1;
MP2_VARopt.quality = 2;
MP2_VARopt.FigSize = [26,12];
MP2_VARopt.firstdate = MP2_datesnum(1);
MP2_VARopt.frequency = 'm';
MP2_VAR.IV = MP2_IV;
MP2_VARopt.ident = 'iv';
MP2_VARopt.snames = {'Odyssean FG Shock','','','','','','','','','','',''};
MP2_VARopt.method = 'wild';
MP2_VARopt.pctg   = 90;

% Compute IRFs
[MP2_IR, MP2_VAR] = VARir(MP2_VAR, MP2_VARopt);
[MP2_IRinf, MP2_IRsup, MP2_IRmed, MP2_IRbar] = VARirband(MP2_VAR, MP2_VARopt);
VARirplot(MP2_IRbar, MP2_VARopt, MP2_IRinf, MP2_IRsup);

% Plot IRFs
FigSize(26, 12);
for ii = 1:MP2_Xnvar
    subplot(4, 3, ii);
    PlotSwathe(MP2_IRbar(:,ii), [MP2_IRinf(:,ii) MP2_IRsup(:,ii)], SwatheOpt); 
    hold on;
    plot(1:MP2_VARopt.nsteps, zeros(MP2_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
    plot(1:MP2_VARopt.nsteps, MP2_IRbar(:,ii), 'Color', cmap(1,:), 'LineWidth', 2);
    title([MP2_Xvnames_long{ii} ' to ' MP2_VARopt.snames{1}], 'FontWeight', 'bold', 'FontSize', 10); 
    xlim([1 MP2_VARopt.nsteps]);
    grid on;
    set(gca, 'Layer', 'top');
end
SaveFigure('figures/cnf_disasters_v10/IRFs_u2', 2);
clf('reset')

%% Counterfactuals
% Regress CPI Food of disasters_IRbar on CPI Food of MP_IRbar
beta = - (MP2_IRbar(:,11) \ data.coeff12(:,3));  % negative to minimize combination
% You get: combined = disasters + beta * MP2 ≈ 0
reg2 = data.coeff12(:,3) + beta * MP2_IRbar(:,11);
plot(reg2)

% Plot IRFs
FigSize(26, 12);
% Assume cmap, SwatheOpt, MP2_IRbar, MP2_IRinf, MP2_IRsup, MP2_VARopt are already defined
% CPI Food
PlotSwathe(data.coeff12(:,3), [data.coeff12_low1(:,3)  data.coeff12_high1(:,3)], SwatheOpt); 
hold on;
plot(1:MP2_VARopt.nsteps, zeros(MP2_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
plot(1:MP2_VARopt.nsteps, data.coeff12(:,3), 'Color', cmap(1,:), 'LineWidth', 2);
plot(1:MP2_VARopt.nsteps, reg2 , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');  % red dashed difference
title([MP2_Xvnames_long{11} ' to disaster shock'], 'FontWeight', 'bold', 'FontSize', 10); 
xlim([1 MP2_VARopt.nsteps]);
SaveFigure('figures/cnf_disasters_v10/cnf_u2', 2);
clf('reset')


%% SECTION 4: MONETARY POLICY (MP) INSTRUMENT VAR FOR LSAP
%------------------------------------------------------------------
[MP3_xlsdata, MP3_xlstext] = xlsread('data/data_JAR_u3.xlsx','Sheet1');
MP3_dates = MP3_xlstext(3:end,1);
MP3_datesnum = Date2Num(MP3_dates, 'm');
MP3_vnames_long = MP3_xlstext(1,2:end);
MP3_vnames = MP3_xlstext(2,2:end);
MP3_nvar = length(MP3_vnames);
MP3_data = Num2NaN(MP3_xlsdata);

% Structure
for ii = 1:MP3_nvar
    MP3_DATA.(MP3_vnames{ii}) = MP3_data(:,ii);
end
MP3_nobs = size(MP3_data,1);

% Log level
log_vars = {'CPIAUCSL','HousePr','StockPr','VIX','INDPRO','CCPI', 'CPIUFDSL'};

for i = 1:length(log_vars)
    MP3_var = log_vars{i};
    MP3_DATA.(MP3_var) = 100 * log(MP3_DATA.(MP3_var));
end

% VAR setup
MP3_Xvnames_long = {'1Y Treasury Bond';'INDPRO';'CPI';'UNRATE';'House Price';'Stock Price';'VIX';'GFC';'EBP';'Commodity Price';'CPI Food'};
MP3_Xvnames = {'DGS1';'INDPRO';'CPIAUCSL';'UNRATE';'HousePr';'StockPr';'VIX';'GFC';'EBP';'CCPI';'CPIUFDSL'};
MP3_Xnvar = length(MP3_Xvnames);
MP3_X = nan(MP3_nobs,MP3_Xnvar);
for ii=1:MP3_Xnvar
    MP3_X(:,ii) = MP3_DATA.(MP3_Xvnames{ii});
end
MP3_IVvnames = {'u3'};
MP3_IV = MP3_DATA.(MP3_IVvnames{1});

% Estimate VAR
MP3_det = 1;
MP3_nlags = 12;
[MP3_VAR, MP3_VARopt] = VARmodel(MP3_X,MP3_nlags,MP3_det);
MP3_VARopt.vnames = MP3_Xvnames_long;
MP3_VARopt.nsteps = 36;
MP3_VARopt.impact    = 1;
MP3_VARopt.quality = 2;
MP3_VARopt.FigSize = [26,12];
MP3_VARopt.firstdate = MP3_datesnum(1);
MP3_VARopt.frequency = 'm';
MP3_VAR.IV = MP3_IV;
MP3_VARopt.ident = 'iv';
MP3_VARopt.ndraws= 1000;
MP3_VARopt.snames = {'LSAP Shock','','','','','','','','','','',''};
MP3_VARopt.method = 'wild';
MP3_VARopt.pctg   = 90;

% Compute IRFs
[MP3_IR, MP3_VAR] = VARir(MP3_VAR, MP3_VARopt);
[MP3_IRinf, MP3_IRsup, MP3_IRmed, MP3_IRbar] = VARirband(MP3_VAR, MP3_VARopt);
VARirplot(MP3_IRbar, MP3_VARopt, MP3_IRinf, MP3_IRsup);

% Plot IRFs
FigSize(26, 12);
for ii = 1:MP3_Xnvar
    subplot(4, 3, ii);
    PlotSwathe(MP3_IRbar(:,ii), [MP3_IRinf(:,ii) MP3_IRsup(:,ii)], SwatheOpt); 
    hold on;
    plot(1:MP3_VARopt.nsteps, zeros(MP3_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
    plot(1:MP3_VARopt.nsteps, MP3_IRbar(:,ii), 'Color', cmap(1,:), 'LineWidth', 2);
    title([MP3_Xvnames_long{ii} ' to ' MP3_VARopt.snames{1}], 'FontWeight', 'bold', 'FontSize', 10); 
    xlim([1 MP3_VARopt.nsteps]);
    grid on;
    set(gca, 'Layer', 'top');
end
SaveFigure('figures/cnf_disasters_v10/IRFs_u3', 2);
clf('reset')

%% Counterfactuals
% Regress CPI Food of disasters_IRbar on CPI Food of MP_IRbar
beta = - (MP3_IRbar(:,11) \ data.coeff12(:,3));  % negative to minimize combination
% You get: combined = disasters + beta * MP ≈ 0
reg3 = data.coeff12(:,3) + beta * MP3_IRbar(:,11);
plot(reg3)

% Plot IRFs
FigSize(26, 12);
% Assume cmap, SwatheOpt, MP_IRbar, MP_IRinf, MP_IRsup, MP_VARopt are already defined
% 1Y Tresury Bond
subplot(2, 3, 1);
PlotSwathe(data.coeff7(:,3), [data.coeff7_low1(:,3)  data.coeff7_high1(:,3)], SwatheOpt); 
hold on;
plot(1:MP3_VARopt.nsteps, zeros(MP3_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
plot(1:MP3_VARopt.nsteps, data.coeff7(:,3), 'Color', cmap(1,:), 'LineWidth', 2);
plot(1:MP3_VARopt.nsteps, reg3 , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');  % red dashed difference
title([MP3_Xvnames_long{1} ' to disaster shock'], 'FontWeight', 'bold', 'FontSize', 10); 
xlim([1 MP3_VARopt.nsteps]);
% INDPRO
subplot(2, 3, 2);
PlotSwathe(data.coeff9(:,3), [data.coeff9_low1(:,3)  data.coeff9_high1(:,3)], SwatheOpt); 
hold on;
plot(1:MP3_VARopt.nsteps, zeros(MP3_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
plot(1:MP3_VARopt.nsteps, data.coeff9(:,3), 'Color', cmap(1,:), 'LineWidth', 2);
plot(1:MP3_VARopt.nsteps, reg3 , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');  % red dashed difference
title([MP3_Xvnames_long{2} ' to disaster shock'], 'FontWeight', 'bold', 'FontSize', 10); 
xlim([1 MP3_VARopt.nsteps]);
% US CPI
subplot(2, 3, 3);
PlotSwathe(data.coeff3(:,3), [data.coeff3_low1(:,3)  data.coeff3_high1(:,3)], SwatheOpt); 
hold on;
plot(1:MP3_VARopt.nsteps, zeros(MP3_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
plot(1:MP3_VARopt.nsteps, data.coeff3(:,3), 'Color', cmap(1,:), 'LineWidth', 2);
plot(1:MP3_VARopt.nsteps, reg3 , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');  % red dashed difference
title([MP3_Xvnames_long{3} ' to disaster shock'], 'FontWeight', 'bold', 'FontSize', 10); 
xlim([1 MP3_VARopt.nsteps]);
% CPI Food
subplot(2, 3, 4);
PlotSwathe(data.coeff12(:,3), [data.coeff12_low1(:,3)  data.coeff12_high1(:,3)], SwatheOpt); 
hold on;
plot(1:MP3_VARopt.nsteps, zeros(MP3_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
plot(1:MP3_VARopt.nsteps, data.coeff12(:,3), 'Color', cmap(1,:), 'LineWidth', 2);
plot(1:MP3_VARopt.nsteps, reg3 , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');  % red dashed difference
title([MP3_Xvnames_long{11} ' to disaster shock'], 'FontWeight', 'bold', 'FontSize', 10); 
xlim([1 MP3_VARopt.nsteps]);
% UNRATE
subplot(2, 3, 5);
PlotSwathe(data.coeff2(:,3), [data.coeff2_low1(:,3)  data.coeff2_high1(:,3)], SwatheOpt); 
hold on;
plot(1:MP3_VARopt.nsteps, zeros(MP3_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
plot(1:MP3_VARopt.nsteps, data.coeff2(:,3), 'Color', cmap(1,:), 'LineWidth', 2);
plot(1:MP3_VARopt.nsteps, reg3 , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');  % red dashed difference
title([MP3_Xvnames_long{4} ' to disaster shock'], 'FontWeight', 'bold', 'FontSize', 10); 
xlim([1 MP3_VARopt.nsteps]);
% House Price
subplot(2, 3, 6);
PlotSwathe(data.coeff4(:,3), [data.coeff4_low1(:,3)  data.coeff4_high1(:,3)], SwatheOpt); 
hold on;
plot(1:MP3_VARopt.nsteps, zeros(MP3_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
plot(1:MP3_VARopt.nsteps, data.coeff4(:,3), 'Color', cmap(1,:), 'LineWidth', 2);
plot(1:MP3_VARopt.nsteps, reg3 , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');  % red dashed difference
title([MP3_Xvnames_long{5} ' to disaster shock'], 'FontWeight', 'bold', 'FontSize', 10); 
xlim([1 MP3_VARopt.nsteps]);
SaveFigure('figures/cnf_disasters_v10/cnf_u3_FINAL', 2);
clf('reset')

%% SECTION 5: MONETARY POLICY (MP) INSTRUMENT FOR Info Effect
%------------------------------------------------------------------
[MP4_xlsdata, MP4_xlstext] = xlsread('data/data_JAR_u4.xlsx','Sheet1');
MP4_dates = MP4_xlstext(3:end,1);
MP4_datesnum = Date2Num(MP4_dates, 'm');
MP4_vnames_long = MP4_xlstext(1,2:end);
MP4_vnames = MP4_xlstext(2,2:end);
MP4_nvar = length(MP4_vnames);
MP4_data = Num2NaN(MP4_xlsdata);

% Structure
for ii = 1:MP4_nvar
    MP4_DATA.(MP4_vnames{ii}) = MP4_data(:,ii);
end
MP4_nobs = size(MP4_data,1);

for i = 1:length(log_vars)
    MP4_var = log_vars{i};
    MP4_DATA.(MP4_var) = 100 * log(MP4_DATA.(MP4_var));
end

% VAR setup
MP4_Xvnames_long = {'1Y Treasury Bond';'INDPRO';'CPI';'UNRATE';'House Price';'Stock Price';'VIX';'GFC';'EBP';'Commodity Price';'CPI Food'};
MP4_Xvnames = {'DGS1';'INDPRO';'CPIAUCSL';'UNRATE';'HousePr';'StockPr';'VIX';'GFC';'EBP';'CCPI';'CPIUFDSL'};
MP4_Xnvar = length(MP4_Xvnames);
MP4_X = nan(MP4_nobs,MP4_Xnvar);
for ii=1:MP4_Xnvar
    MP4_X(:,ii) = MP4_DATA.(MP4_Xvnames{ii});
end
MP4_IVvnames = {'u4'};
MP4_IV = MP4_DATA.(MP4_IVvnames{1});

% Estimate VAR
MP4_det = 1;
MP4_nlags = 12;
[MP4_VAR, MP4_VARopt] = VARmodel(MP4_X,MP4_nlags,MP4_det);
MP4_VARopt.vnames = MP4_Xvnames_long;
MP4_VARopt.nsteps = 36;
MP4_VARopt.impact    = 1;
MP4_VARopt.quality = 2;
MP4_VARopt.FigSize = [26,12];
MP4_VARopt.firstdate = MP4_datesnum(1);
MP4_VARopt.frequency = 'm';
MP4_VAR.IV = MP4_IV;
MP4_VARopt.ident = 'iv';
MP4_VARopt.snames = {'Delphic shock','','','','','','','','','','',''};
MP4_VARopt.method = 'wild';
MP4_VARopt.pctg   = 90;

% Compute IRFs
[MP4_IR, MP4_VAR] = VARir(MP4_VAR, MP4_VARopt);
[MP4_IRinf, MP4_IRsup, MP4_IRmed, MP4_IRbar] = VARirband(MP4_VAR, MP4_VARopt);
VARirplot(MP4_IRbar, MP4_VARopt, MP4_IRinf, MP4_IRsup);

% Plot IRFs
cmap = lines;
FigSize(26, 12);
SwatheOpt = PlotSwatheOption;
SwatheOpt.swathecol = [0.85 0.85 0.85];
SwatheOpt.linecol = [0.2 0.2 0.2];
for ii = 1:MP4_Xnvar
    subplot(4, 3, ii);
    PlotSwathe(MP4_IRbar(:,ii), [MP4_IRinf(:,ii) MP4_IRsup(:,ii)], SwatheOpt); 
    hold on;
    plot(1:MP4_VARopt.nsteps, zeros(MP4_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
    plot(1:MP4_VARopt.nsteps, MP4_IRbar(:,ii), 'Color', cmap(1,:), 'LineWidth', 2);
    title([MP4_Xvnames_long{ii} ' to ' MP4_VARopt.snames{1}], 'FontWeight', 'bold', 'FontSize', 10); 
    xlim([1 MP4_VARopt.nsteps]);
    grid on;
    set(gca, 'Layer', 'top');
end
SaveFigure('figures/cnf_disasters_v10/IRFs_u4', 2);
clf('reset')

%% Counterfactuals
% Regress CPI Food of disasters_IRbar on CPI Food of MP_IRbar
beta = - (MP4_IRbar(:,11) \ data.coeff12(:,3));  % negative to minimize combination
% You get: combined = disasters + beta * MP ≈ 0
reg4 = data.coeff12(:,3) + beta * MP4_IRbar(:,11);
plot(reg4)

% Plot IRFs
FigSize(26, 12);
% Assume cmap, SwatheOpt, MP_IRbar, MP_IRinf, MP_IRsup, MP_VARopt are already defined
% CPI Food
PlotSwathe(data.coeff12(:,3), [data.coeff12_low1(:,3)  data.coeff12_high1(:,3)], SwatheOpt); 
hold on;
plot(1:MP4_VARopt.nsteps, zeros(MP4_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
plot(1:MP4_VARopt.nsteps, data.coeff12(:,3), 'Color', cmap(1,:), 'LineWidth', 2);
plot(1:MP4_VARopt.nsteps, reg4 , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');  % red dashed difference
title([MP4_Xvnames_long{11} ' to disaster shock'], 'FontWeight', 'bold', 'FontSize', 10); 
xlim([1 MP4_VARopt.nsteps]);
SaveFigure('figures/cnf_disasters_v10/cnf_u4', 2);
clf('reset')

%% SECTION 6: Counterfactuals
%------------------------------------------------------------------
% Minimising the effect of MP shock 
beta = - (MP_IRbar(:,11) \ data.coeff12(:,3));  % negative to minimize combination
% You get: combined = disasters + beta * MP ≈ 0
reg_u1 = data.coeff12(:,3) + beta * MP_IRbar(:,11);

beta = - (MP2_IRbar(:,11) \ data.coeff12(:,3));  % negative to minimize combination
% You get: combined = disasters + beta * MP ≈ 0
reg_u2 = data.coeff12(:,3) + beta * MP2_IRbar(:,11);

beta = - (MP3_IRbar(:,11) \ data.coeff12(:,3));  % negative to minimize combination
% You get: combined = disasters + beta * MP ≈ 0
reg_u3 = data.coeff12(:,3) + beta * MP3_IRbar(:,11);

beta = - (MP4_IRbar(:,11) \ data.coeff12(:,3));  % negative to minimize combination
% You get: combined = disasters + beta * MP ≈ 0
reg_u4 = data.coeff12(:,3) + beta * MP4_IRbar(:,11);

% Plot regs
FigSize(26, 12);
plot(1:MP_VARopt.nsteps, reg_u1 , 'LineWidth', 2); 
hold on;
plot(1:MP_VARopt.nsteps, reg_u2 ,'LineWidth', 2);  
plot(1:MP_VARopt.nsteps, reg_u3 , 'LineWidth', 2);  
plot(1:MP_VARopt.nsteps, reg_u4 , 'LineWidth', 2);  
legend('Standard shock (u1)', 'Odyssean guidance (u2)', 'LSAP shock (u3)', 'Delphic guidance (u4)');
xlim([1 MP_VARopt.nsteps]);
SaveFigure('figures/cnf_disasters_v10/regs', 2);
clf('reset')

% Counterfactuals
X = [MP3_IRbar(:,11), MP4_IRbar(:,11)];
% Target (negated disasters)
y = data.coeff12(:,3);
% Solve regression
beta2 = -(X \ y);  % returns [beta1; beta2]

% Construct the linear combination
reg13 = data.coeff12(:,3) + beta2(1) * MP_IRbar(:,11) + beta2(2) * MP3_IRbar(:,11);
plot(reg13)
title('Minimized Linear Combination (3 Series)')
xlabel('Time')
ylabel('Combined Series')
clf('reset')

% Plot IRFs
FigSize(26, 12);
% Assume cmap, SwatheOpt, MP_IRbar, MP_IRinf, MP_IRsup, MP_VARopt are already defined
% CPI Food
PlotSwathe(data.coeff12(:,3), [data.coeff12_low1(:,3)  data.coeff12_high1(:,3)], SwatheOpt); 
hold on;
plot(1:MP_VARopt.nsteps, zeros(MP_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
plot(1:MP_VARopt.nsteps, data.coeff12(:,3), 'Color', cmap(1,:), 'LineWidth', 2);
plot(1:MP_VARopt.nsteps, reg13 , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');  % red dashed difference
title([MP_Xvnames_long{11} ' to disaster shock'], 'FontWeight', 'bold', 'FontSize', 10); 
xlim([1 MP_VARopt.nsteps]);
SaveFigure('figures/cnf_disasters_v10/cnf_u1u3', 2);
clf('reset')
