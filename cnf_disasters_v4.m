% this code modifies v3 in following ways:
% - change VAR for LP for natural disaster shock 

%% 0. PRELIMINARIES
%------------------------------------------------------------------
clear all; close all; clc
warning off all
format short g

%% Add path to toolbox

% Create output folders
mkdir('figures/cnf_disasters_v4')

% commiting a change

%% SECTION 1: disasters INSTRUMENT VAR
%------------------------------------------------------------------
% Load data
[disasters_xlsdata, disasters_xlstext] = xlsread('data/data_disasters_v1.xlsx','Sheet1');
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
    title(disasters_vnames_long(ii)); 
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
log_vars = {'CPIAUCSL','HousePr','StockPr','VIX','INDPRO','CCPI'};

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
disaster = data.STORMSFLOODS;

time = datetime(2000,1,1):calmonths(1):datetime(2019,12,1);

figure;
bar(time, disaster)
SaveFigure('figures/cnf_disasters_v4/EMDAT',2)

% Disaster shock: contemporaneous + 9 lags 
numLags=9;
disaster_lags = lagmatrix(disaster, 1:numLags);
disaster_shock = [disaster, disaster_lags];
disaster_shock = disaster_shock(numLags+1:end, :);

% Lagged control variables
numLagsCV=3;
% List of control variables (FD)
lagged_vars = {'UNRATE_FD','CPIAUCSL_FD','HousePr_FD','StockPr_FD','VIX_FD',...
                'DGS1_FD','EBP_FD','INDPRO_FD','GFC','CCPI_FD'};

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
response_vars = {'UNRATE','CPIAUCSL','HousePr','StockPr','VIX','DGS1','GFC','INDPRO','EBP','CCPI'};

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

horizons=36; 

%% Long-run differences
dependent_indices = 2:11;  % Corresponds to Y2 to Y11
confidence_levels = [1, 1.645];  % 68% and 90% confidence intervals

for idx = dependent_indices
    Y = data.(['Y' num2str(idx)]);

    % Select appropriate X matrix
    if ismember(idx, [9, 10, 11])
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
    struct('title', 'UNRATE', 'base', 'coeff2'), ...    
    struct('title', 'CPIAUCSL', 'base', 'coeff3'), ...
    struct('title', 'House Price', 'base', 'coeff4'), ...
    struct('title', 'Stock Price', 'base', 'coeff5'), ...
    struct('title', 'VIX', 'base', 'coeff6'), ...    
    struct('title', 'DGS1', 'base', 'coeff7'), ...    
    struct('title', 'INDPRO', 'base', 'coeff9'), ...
    struct('title', 'EBP', 'base', 'coeff10'), ...
    struct('title', 'CCPI', 'base', 'coeff11')
};

for i = 1:length(plot_configs)
    subplot(3, 3, i);
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
SaveFigure('figures/cnf_disasters_v4/IRFs_disasters', 2);
clf('reset')

%% SECTION 2: MONETARY POLICY (MP) INSTRUMENT VAR
%------------------------------------------------------------------
[MP_xlsdata, MP_xlstext] = xlsread('data/data_MAR_d.xlsx','Sheet1');
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

% VAR setup
% MAR21 use VAR(12) and even with 3 variables VAR (IP, CPI, GS1) got significant results
% EBP is in log in the paper!
% puzzling is that commodity (energy here) price index is not negative 
MP_Xvnames_long = {'1Y Treasury Bond';'INDPRO';'CPI';'UNRATE';'House Price';'Stock Price';'VIX';'GFC';'EBP';'Commodity Price'};
MP_Xvnames = {'DGS1';'INDPRO';'CPIAUCSL';'UNRATE';'HousePr';'StockPr';'VIX';'GFC';'EBP';'CCPI'};
MP_Xnvar = length(MP_Xvnames);
MP_X = nan(MP_nobs,MP_Xnvar);
for ii=1:MP_Xnvar
    MP_X(:,ii) = MP_DATA.(MP_Xvnames{ii});
end
MP_IVvnames = {'MAR'};
MP_IV = MP_DATA.(MP_IVvnames{1});

% Estimate VAR
MP_det = 1;
MP_nlags = 12;
[MP_VAR, MP_VARopt] = VARmodel(MP_X,MP_nlags,MP_det);
MP_VARopt.vnames = MP_Xvnames_long;
MP_VARopt.nsteps = 36;
%MP_VARopt.impact    = 1;
MP_VARopt.quality = 2;
MP_VARopt.FigSize = [26,12];
MP_VARopt.firstdate = MP_datesnum(1);
MP_VARopt.frequency = 'm';
MP_VAR.IV = MP_IV;
MP_VARopt.ident = 'iv';
MP_VARopt.snames = {'MP Shock','','','','','','','','','',''};
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
SaveFigure('figures/cnf_disasters_v4/IRFs_MAR', 2);
clf('reset')

%% Counterfactuals
% Regress 1Y Tresury Bond of disasters_IRbar on 1Y Tresury Bond of MP_IRbar
beta = - (MP_IRbar(:,1) \ data.coeff7(:,3));  % negative to minimize combination
% You get: combined = disasters + beta * MP ≈ 0
reg1 = data.coeff7(:,3) + beta * MP_IRbar(:,1);

% Plot IRFs
FigSize(26, 12);
% Assume cmap, SwatheOpt, MP_IRbar, MP_IRinf, MP_IRsup, MP_VARopt are already defined
% 1Y Tresury Bond
PlotSwathe(data.coeff7(:,3), [data.coeff7_low1(:,3)  data.coeff7_high1(:,3)], SwatheOpt); 
hold on;
plot(1:MP_VARopt.nsteps, zeros(MP_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
plot(1:MP_VARopt.nsteps, data.coeff7(:,3), 'Color', cmap(1,:), 'LineWidth', 2);
plot(1:MP_VARopt.nsteps, reg1 , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');  % red dashed difference
title([MP_Xvnames_long{1} ' to disaster shock'], 'FontWeight', 'bold', 'FontSize', 10); 
xlim([1 MP_VARopt.nsteps]);
SaveFigure('figures/cnf_disasters_v4/cnf_MAR', 2);
clf('reset')

%% SECTION 3: MONETARY POLICY (MP) INSTRUMENT VAR
%------------------------------------------------------------------
[MP2_xlsdata, MP2_xlstext] = xlsread('data/data_BRW_d.xlsx','Sheet1');
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

% VAR setup
MP2_Xvnames_long = {'1Y Treasury Bond';'INDPRO';'CPI';'UNRATE';'House Price';'Stock Price';'VIX';'GFC';'EBP';'Commodity Price'};
MP2_Xvnames = {'DGS1';'INDPRO';'CPIAUCSL';'UNRATE';'HousePr';'StockPr';'VIX';'GFC';'EBP';'CCPI'};
MP2_Xnvar = length(MP2_Xvnames);
MP2_X = nan(MP2_nobs,MP2_Xnvar);
for ii=1:MP2_Xnvar
    MP2_X(:,ii) = MP2_DATA.(MP2_Xvnames{ii});
end
MP2_IVvnames = {'BRW'};
MP2_IV = MP2_DATA.(MP2_IVvnames{1});

% Estimate VAR
MP2_det = 1;
MP2_nlags = 12;
[MP2_VAR, MP2_VARopt] = VARmodel(MP2_X,MP2_nlags,MP2_det);
MP2_VARopt.vnames = MP2_Xvnames_long;
MP2_VARopt.nsteps = 36;
MP2_VARopt.quality = 2;
MP2_VARopt.FigSize = [26,12];
MP2_VARopt.firstdate = MP2_datesnum(1);
MP2_VARopt.frequency = 'm';
MP2_VAR.IV = MP2_IV;
MP2_VARopt.ident = 'iv';
MP2_VARopt.ndraws= 1000;
MP2_VARopt.snames = {'MP2 Shock','','','','','','','','','',''};
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
SaveFigure('figures/cnf_disasters_v4/IRFs_BRW', 2);
clf('reset')

%% Counterfactuals
% Regress 1Y Tresury Bond of disasters_IRbar on 1Y Tresury Bond of MP_IRbar
beta = - (MP2_IRbar(:,1) \ data.coeff7(:,3));  % negative to minimize combination
% You get: combined = disasters + beta * MP ≈ 0
reg2 = data.coeff7(:,3) + beta * MP2_IRbar(:,1);

% Plot IRFs
FigSize(26, 12);
% Assume cmap, SwatheOpt, MP2_IRbar, MP2_IRinf, MP2_IRsup, MP2_VARopt are already defined
% 1Y Tresury Bond
PlotSwathe(data.coeff7(:,3), [data.coeff7_low1(:,3)  data.coeff7_high1(:,3)], SwatheOpt); 
hold on;
plot(1:MP2_VARopt.nsteps, zeros(MP2_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
plot(1:MP2_VARopt.nsteps, data.coeff7(:,3), 'Color', cmap(1,:), 'LineWidth', 2);
plot(1:MP2_VARopt.nsteps, reg2 , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');  % red dashed difference
title([MP2_Xvnames_long{1} ' to disaster shock'], 'FontWeight', 'bold', 'FontSize', 10); 
xlim([1 MP2_VARopt.nsteps]);
SaveFigure('figures/cnf_disasters_v4/cnf_BRW', 2);
clf('reset')

% Plot regs
FigSize(26, 12);
plot(1:MP_VARopt.nsteps, reg1 , 'LineWidth', 2); 
hold on;
plot(1:MP_VARopt.nsteps, reg2 ,'LineWidth', 2);   
legend('MAR', 'BRW');
xlim([1 MP_VARopt.nsteps]);
SaveFigure('figures/cnf_disasters_v4/regs', 2);

% Linear combination
X = [MP_IRbar(:,1), MP2_IRbar(:,1)];
% Target (negated disasters)
y = data.coeff7(:,3);
% Solve regression
beta2 = -(X \ y);  % returns [beta1; beta2]

% Construct the linear combination
reg3 = data.coeff7(:,3) + beta2(1) * MP_IRbar(:,1) + beta2(2) * MP2_IRbar(:,1);
title('Minimized Linear Combination (3 Series)')
xlabel('Time')
ylabel('Combined Series')

% Plot IRFs
FigSize(26, 12);
% Assume cmap, SwatheOpt, MP_IRbar, MP_IRinf, MP_IRsup, MP_VARopt are already defined
% 1Y Tresury Bond
PlotSwathe(data.coeff7(:,3), [data.coeff7_low1(:,3)  data.coeff7_high1(:,3)], SwatheOpt); 
hold on;
plot(1:MP_VARopt.nsteps, zeros(MP_VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
plot(1:MP_VARopt.nsteps, data.coeff7(:,3), 'Color', cmap(1,:), 'LineWidth', 2);
plot(1:MP_VARopt.nsteps, reg3 , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');  % red dashed difference
title([MP_Xvnames_long{1} ' to disaster shock'], 'FontWeight', 'bold', 'FontSize', 10); 
xlim([1 MP_VARopt.nsteps]);
SaveFigure('figures/cnf_disasters_v4/cnf_2MP', 2);
clf('reset')