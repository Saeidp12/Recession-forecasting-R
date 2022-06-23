%% mainRecForecast.m
% DESCRIPTION: 
% Main entry point for recession forecasting framework.
%
% AUTHOR: 
% Aaron Smalter Hall, April 2017
% Federal Reserve Bank of Kansas City
%
% RELATED CITATION:
% Recession Forecasting Using Bayesian Classification
% Troy Davig and Aaron Smalter Hall
% August 2016; Revised February 2017
% RWP 16-06
% https://dx.doi.org/10.18651/RWP2016-06

%% Clear all existing variables
clear;

%% Define settings for forecast

% -------------------------------
% Data file and variable settings
% -------------------------------
rec_file = 'NBER_recessions.csv';
macro_file = 'macro_data.csv';
selected_macro_vars = {'NAPMPI','PAYEMS_RT','S&P 500','T10YFFM'};

% -------------------------------
% Lag structure settings
% -------------------------------
n_lag = 10;

% -------------------------------
% Forecast settings
% -------------------------------
forecast_horizons = [0,3,6,9,12];
n_horz = length(forecast_horizons);

forecast_blind = 18;
use_alt_labels = true; % use "recession within horizon" labels instead of "recession at exactly horizon"

forecast_start_date = 1973;
forecast_end_date = Inf; % using 'Inf' for end date will stop at the last date we have data for

% -------------------------------
% Model settings
% -------------------------------
% use_model = 1; % use logistic regression model
% use_model = 2; % use naive bayes model
use_model = 3; % use modified naive bayes model

suppress_warns = true;

% -------------------------------
% Performance evaluation settings
% -------------------------------
use_timedep_weights = true; % use timing-dependent error weights instead of uniform weights

%% Read NBER recession data from file
rec_data_raw = readtable(rec_file);
rec_labels = rec_data_raw.Recession;
[rec_years,rec_months] = datevec(rec_data_raw.Date);
rec_dates = rec_years+(rec_months-1)/12;

%% Read macroeconomic data from file
generate_freddata;
macro_data = data;
macro_dates = dates;
macro_vars = names;

%% Define macro variables to be used for forecasting
selected_var_idx = find(ismember(macro_vars,selected_macro_vars));
macro_data_selected = macro_data(:,selected_var_idx);

%% Define and generate lags of macro variables
macro_data_lagged = nan(size(macro_data_selected,1),(n_lag+1)*size(macro_data_selected,2));
for i=1+n_lag:size(macro_data_selected,1)
    data_lagged = macro_data_selected(i-n_lag:i,:);
    macro_data_lagged(i,:) = data_lagged(:);
end

%% Define final macro data and dates
macro_data_final = macro_data_lagged(n_lag+1:end,:);
macro_dates_final = macro_dates(n_lag+1:end,:);

%% Find common date range for recession and macro data
% WARNING: Dates in both recession and macro data must be consecutive
% and without gaps, there are no checks for catching or recovering in case
% date ranges are ill-formed.

common_dates = intersect(single(macro_dates_final),single(rec_dates));
macro_common_idx = find(ismember(single(macro_dates_final),common_dates));
rec_common_idx = find(ismember(single(rec_dates),common_dates));

macro_data_common = macro_data_final(macro_common_idx,:);
macro_dates_common = macro_dates_final(macro_common_idx,:);

rec_labels_common = rec_labels(rec_common_idx,:);
rec_dates_common = rec_dates(rec_common_idx,:);

[n_obs,n_vars] = size(macro_data_common);
        % test_y = rec_labels_common(test_idx,:);
%% Create data structures and indices for forecast
forecast_end_date = min(forecast_end_date,macro_dates_common(end,:));
forecast_dates = (macro_dates_common(1,:):1/12:macro_dates_common(end,:)+max(forecast_horizons)/12)';
forecast_start_idx = find(single(forecast_start_date)==single(macro_dates_common));
forecast_end_idx = find(single(forecast_end_date)==single(macro_dates_common));
forecast_probs = nan(n_obs+max(forecast_horizons),n_horz);

%% Set warning suppression
if suppress_warns
    warning('off','stats:glmfit:PerfectSeparation');
    warning('off','stats:glmfit:IterationLimit');
else
    warning('on','stats:glmfit:PerfectSeparation');
    warning('on','stats:glmfit:IterationLimit');
end

%% Iterate through forecast dates and horizons, generating forecast probabilities along the way
disp('---------------------------');
disp('Starting Recession Forecast');
disp('---------------------------');
tic
for forecast_i=forecast_start_idx:forecast_end_idx
    pct_progress = 100*(forecast_i-forecast_start_idx)/(forecast_end_idx-forecast_start_idx);
    if mod(100*(forecast_i-forecast_start_idx)/round((forecast_end_idx-forecast_start_idx)/10),10)==0
        fprintf('Forecasting... %d/%d \t (%.2f%%)  [%.2fs]\n',forecast_i-forecast_start_idx,forecast_end_idx-forecast_start_idx,pct_progress,toc);
    end
    
    for horizon_i=1:n_horz
        forecast_horizon = forecast_horizons(horizon_i);
        
        train_idx = 1:forecast_i-1-forecast_blind-forecast_horizon;
        train_X = macro_data_common(train_idx,:);
        train_y = rec_labels_common(train_idx+forecast_horizon,:);
        
        test_idx = forecast_i;
        test_X = macro_data_common(test_idx,:);
        
        if use_alt_labels
            alt_train_y = train_y;
            for d_i = (1+forecast_horizon):length(train_y)
                if train_y(d_i-forecast_horizon) == 0 && sum(train_y((d_i-forecast_horizon):d_i))
                    alt_train_y(d_i) = 1;
                end
            end
            train_y = alt_train_y;
        end
        
        if use_model == 1 % use logistic regression model
            model = fitglm(train_X,train_y,'Distribution','binomial','Link','logit');
            prob = predict(model,test_X);
            probs = [1-prob,prob];
        elseif use_model == 2 % use naive bayes model
            model = fitcnb(train_X,train_y,'ClassNames',[0;1]);
            [~,probs] = predict(model,test_X);
        else % default to MNB if no model selected
            last_pos_prob = forecast_probs(forecast_i+forecast_horizon-1,horizon_i);
            last_neg_prob = 1-last_pos_prob;

            pos_after_neg_prob = sum(train_y(find(train_y(1:end-1)==0)+1)==1)/sum(train_y(1:end-1)==0);
            neg_after_neg_prob = sum(train_y(find(train_y(1:end-1)==0)+1)==0)/sum(train_y(1:end-1)==0);

            pos_after_pos_prob = sum(train_y(find(train_y(1:end-1)==1)+1)==1)/sum(train_y(1:end-1)==1);
            neg_after_pos_prob = sum(train_y(find(train_y(1:end-1)==1)+1)==0)/sum(train_y(1:end-1)==1);

            pos_prior = pos_after_neg_prob*last_neg_prob+pos_after_pos_prob*last_pos_prob;
            neg_prior = neg_after_neg_prob*last_neg_prob+neg_after_pos_prob*last_pos_prob;
            
            modified_priors = [neg_prior,pos_prior];
            if isnan(modified_priors)
                modified_priors = 'empirical';
            end
            
            model = fitcnb(train_X,train_y,'Prior',modified_priors,'ClassNames',[0;1]);
            [~,probs] = predict(model,test_X);
        end
        
        forecast_probs;
    end
end
fprintf('Forecast done! %d/%d \t (%.2f%%)  [%.2fs]\n',forecast_i-forecast_start_idx,forecast_end_idx-forecast_start_idx,pct_progress,toc);
disp('---------------------------');
toc
%% Construct forecast weights for performance calculations
forecast_weights = ones(size(forecast_probs,1),1);
if use_timedep_weights
    label_edges = (rec_labels_common(1:end-1)==0 & rec_labels_common(2:end)==1) | (rec_labels_common(1:end-1)==1 & rec_labels_common(2:end)==0);
    label_edge_idxs = [0;find(label_edges);length(label_edges)+1];

    for i=1:length(label_edge_idxs)-1
        idx_range = label_edge_idxs(i)+1:label_edge_idxs(i+1);
        range_size = length(idx_range);
        weights = normpdf(0:1/(range_size-1):1,0.5,1);
        scaled_weights = (weights-min(weights))./(max(weights)-min(weights));

        if i==1
            scaled_weights(1:round(length(scaled_weights)/2))=max(scaled_weights);
        elseif i==length(label_edge_idxs)-1
            scaled_weights(round(length(scaled_weights)/2):end)=max(scaled_weights);
        end

        forecast_weights(idx_range) = scaled_weights.*range_size./sum(scaled_weights);
    end
end

%% Assemble forecasts into a table for output and review
forecast_years = floor(forecast_dates);
forecast_months = round(mod(forecast_dates,1)*12+1);
forecast_table = table(forecast_years,forecast_months,'VariableNames',{'date_year','date_month'});
for horizon_i=1:n_horz
    forecast_table.(['H' num2str(forecast_horizons(horizon_i))]) = forecast_probs(:,horizon_i);
end
forecast_table.recession = [rec_labels_common;nan(max(forecast_horizons),1)];

fprintf('\nRecent Forecasts:\n');
disp(forecast_table(end-max(forecast_horizons):end,:));
fprintf('*** Forecasts saved to variable: forecast_table\n');

writetable(forecast_table);

%% Calculate and print summary statistics for forecast performance
performance_table = table({'MAE';'MSE';'AUC';'ACC';'F1'},'VariableNames',{'statistic'});
for horizon_i=1:n_horz
    horizon_probs = forecast_probs(:,horizon_i);
    horizon_idx = find(~isnan(horizon_probs));
    
    rec_outcomes = forecast_table.recession;
    rec_idx = find(~isnan(rec_outcomes));
    
    common_idx = intersect(horizon_idx,rec_idx);
    horizon_probs_common = horizon_probs(common_idx,:);
    rec_outcomes_common = rec_outcomes(common_idx,:);
    forecast_weights_common = forecast_weights(common_idx,:);
    
    MAE = mean(abs(horizon_probs_common-rec_outcomes_common).*forecast_weights_common);
    MSE = mean((horizon_probs_common-rec_outcomes_common).^2.*forecast_weights_common);
    
    [~,~,~,AUC] = perfcurve(rec_outcomes_common,horizon_probs_common,1,'Weights',forecast_weights_common);
    
    TP = sum((rec_outcomes_common==1 & horizon_probs_common>0.5).*forecast_weights_common);
    TN = sum((rec_outcomes_common==0 & horizon_probs_common<0.5).*forecast_weights_common);
    FP = sum((rec_outcomes_common==0 & horizon_probs_common>0.5).*forecast_weights_common);
    FN = sum((rec_outcomes_common==1 & horizon_probs_common<0.5).*forecast_weights_common);

    precision = TP/(TP+FP);
    recall = TP/(TP+FN);
    
    F1 = max(0,2*(precision*recall)/(precision+recall));
    ACC = (TP+TN)/(TP+TN+FP+FN);
    
    performance_table.(['H' num2str(forecast_horizons(horizon_i))]) = [MAE;MSE;AUC;ACC;F1];
end
performance_table.AVG = mean(performance_table{:,2:n_horz+1},2);
writetable(performance_table);

fprintf('\nPerformance Table:\n');
disp(performance_table);
fprintf('*** Performance saved to variable: performance_table\n');

%%