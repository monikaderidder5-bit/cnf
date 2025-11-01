% this code modifies v14 in following ways:
% VAR instead of LP

%% 0. PRELIMINARIES
%------------------------------------------------------------------
clear all; close all; clc
warning off all
format short g

%% Add path to toolbox
% addpath('functions')

% Create output folders
mkdir('figures/cnf_disasters_v14')

%% SECTION 1: disasters INSTRUMENT VAR
%------------------------------------------------------------------
% Load data
[xlsdata, xlstext] = xlsread('data/data_affected.xlsx','Sheet1');
dates = xlstext(3:end,1);
datesnum = Date2Num(dates, 'm');
vnames_long = xlstext(1,2:end);
vnames = xlstext(2,2:end);
nvar = length(vnames);
data = Num2NaN(xlsdata);

% Structure
for ii = 1:nvar
    DATA.(vnames{ii}) = data(:,ii);
end
nobs = size(data,1);

% Log level
log_vars = {'CPIAUCSL','HousePr','StockPr','VIX','INDPRO','CCPI', 'CPIUFDSL'};

for i = 1:length(log_vars)
    var = log_vars{i};
    DATA.(var) = 100 * log(DATA.(var));
end

%% Define shock variables and their labels
shock_vars = {'TotalaffectedFlood', 'TotalaffectedStorm', 'TotalaffectedWildfire'};
shock_labels = {'Flood', 'Storm', 'Wildfire'};

for s = 1:length(shock_vars)
    
    shock_var = shock_vars{s};
    shock_label = shock_labels{s};

    %----------------------------------------------
    % Normalize disaster variable (mean of non-zero = 1)
    raw_disaster = DATA.(shock_var);  % Get raw data
    weights_temp = zeros(size(raw_disaster));
    weights = zeros(size(raw_disaster));
    nonzero_idx = raw_disaster > 0;

    if any(nonzero_idx)
        minimum = min(raw_disaster(nonzero_idx));
        avg = mean(raw_disaster(nonzero_idx));

        for i = 1:length(raw_disaster)
            weights_temp(i) = raw_disaster(i) * minimum;
        end

        for i = 1:length(raw_disaster)
            if raw_disaster(i) ~= 0
                weights(i) = (nnz(raw_disaster) / sum(weights_temp)) * weights_temp(i);
            end
        end
    else
        warning('%s contains no non-zero values — skipping normalization.', shock_var);
        weights = raw_disaster;
        avg = 0;
    end

    % Use normalized disaster as new shock variable
    DATA.disaster = weights;
    clean_shock_name = 'disaster';  % name for normalized variable in VAR

    %% SECTION 2: VAR ESTIMATION
    %----------------------------------------------
    Xvnames = {clean_shock_name, 'INDPRO', 'CPIAUCSL', 'DGS1'};
    Xvnames_long = {shock_label, 'INDPRO', 'CPI', 'DGS1'};
    Xnvar = length(Xvnames);
    
    % Create matrix of VAR variables
    X = nan(nobs, Xnvar);
    for ii = 1:Xnvar
        X(:, ii) = DATA.(Xvnames{ii});
    end

    % Make a common sample by removing NaNs
    [X, fo, lo] = CommonSample(X);
    
    % Set VAR parameters
    det = 1;
    nlags = 12;
    
    % Estimate VAR
    [VAR, VARopt] = VARmodel(X, nlags, det);
    VARopt.vnames = Xvnames_long;
    
    % Print VAR results
    [TABLE, beta] = VARprint(VAR, VARopt, 2);

    %% SECTION 3: IDENTIFICATION AND IRFs
    %----------------------------------------------
    VARopt.ident = 'short';
    VARopt.ndraws = 1000;
    VARopt.vnames = Xvnames_long;
    VARopt.nsteps = 48;
    VARopt.impact = 1;
    VARopt.FigSize = [26, 12];
    VARopt.firstdate = datesnum(1);
    VARopt.frequency = 'm';
    VARopt.snames = {[lower(shock_label) ' shock'], '', '', '', ''};
    VARopt.method = 'bs';
    VARopt.pctg = 90;

    % IRF computation
    [IR, VAR] = VARir(VAR, VARopt);
    eps_short = (VAR.B \ VAR.resid')';
    disp(['Correlation matrix of structural shocks for: ' shock_label]);
    disp(corr(eps_short));

    % IR bands and plots
    [IRinf, IRsup, IRmed, IRbar] = VARirband(VAR, VARopt);
    VARirplot(IRbar, VARopt, IRinf, IRsup);

    %% Plot IRFs (custom format)
    cmap = lines;
    FigSize(26, 12);
    SwatheOpt = PlotSwatheOption;
    SwatheOpt.swathecol = [0.85 0.85 0.85];
    SwatheOpt.linecol = [0.2 0.2 0.2];

    sgtitle(sprintf('Impulse: %s \nAvg. of non-zero event: %dK people', ...
        shock_label, round(avg / 1000)));

    for i = 1:Xnvar
        subplot(2,2,i);
        PlotSwathe(IRbar(:,i), [IRinf(:,i), IRsup(:,i)], SwatheOpt); 
        hold on;
        plot(1:VARopt.nsteps, zeros(VARopt.nsteps,1), '--k', 'LineWidth', 0.5);
        plot(1:VARopt.nsteps, IRbar(:,i), 'Color', cmap(1,:), 'LineWidth', 2);
        title([Xvnames_long{i} ' to ' VARopt.snames{1}], 'FontWeight', 'bold', 'FontSize', 10); 
        xlim([1 VARopt.nsteps]);
        grid on;
        set(gca, 'Layer', 'top');
    end

    % Save the figure
    SaveFigure(['figures/cnf_disasters_v14/normalized_IRFs_VAR_affected_' lower(shock_label)], 2);
    clf('reset');  % Reset figure for next loop

end
