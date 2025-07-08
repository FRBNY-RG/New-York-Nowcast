%%% Nowcasting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script produces a nowcast of US real GDP growth for 2023:Q4 
% using the estimated parameters from a dynamic factor model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear workspace and set paths
close all; clear; clc;

addpath('functions')
addpath(['functions' filesep 'general'])
addpath(['functions' filesep 'model'])

% Create function handle for data filenames and directories
fname = @(type, d) [type '_' datestr(d, 'yyyy_mm_dd')];

%% User inputs
date_nowcast_new  = datetime(2023, 09, 29); % this Friday's date
date_estimate_new = datetime(2023, 09, 20); % estimate file for current week

date_nowcast_old  = datetime(2023, 09, 22); % last Friday's date
date_estimate_old = datetime(2023, 09, 20); % estimate file for previous week

date_forecast_old = datetime(2023, 12, 01);
date_forecast_new = datetime(2023, 12, 01);

% Set seed for simulation-based calculations
rng(321)

%% Load model specification and data
spec = load_spec('model_spec_FRED.csv', 1);

Data_nowcast_old = load(['data\' fname('Data', date_nowcast_old)]);
Data_nowcast_new = load(['data\' fname('Data', date_nowcast_new)]);

% Upload initial values for parameters
load('initval.mat'); 

% Upload model settings
settings = load_settings();

% Pair data and specification (nowcasting) 
n = size(spec.SeriesID, 1);

% Find index for GDP (variable to nowcast)
i_now = find( strcmpi('GDPC1',spec.SeriesID));

% Extend datasets to cover forecasting period
Data_nowcast_old.timekey = datetime(Data_nowcast_old.timekey);
date_tmp                 = ((Data_nowcast_old.timekey(end)+calmonths(1)):calmonths:date_forecast_old)';
Data_nowcast_old.timekey = [Data_nowcast_old.timekey; date_tmp];
Data_nowcast_old.data    = [Data_nowcast_old.data; NaN(length(date_tmp), n)];
Data_nowcast_new.timekey = datetime(Data_nowcast_new.timekey);
date_tmp                 = ((Data_nowcast_new.timekey(end)+calmonths(1)):calmonths:date_forecast_new)';
Data_nowcast_new.timekey = [Data_nowcast_new.timekey; date_tmp];
Data_nowcast_new.data    = [Data_nowcast_new.data; NaN(length(date_tmp), n)];

% Extend initial values to match data
T = size(Data_nowcast_new.data, 1);
excess = T -size(initval.sigma, 2);
if excess > 0
    for i = 1:excess
        initval.sigma = [initval.sigma initval.sigma(:, end)];
        initval.s = [initval.s initval.s(:, end)];
    end
end 


%% Nowcasting

fprintf('NOWCAST\n\n')

% Upload estimation output
estimates = load(fname('Estimates', date_estimate_new));
estimates_old = load(fname('Estimates', date_estimate_old));

% Fix location and scale of the data
Y_location = estimates.Y_location;
Y_scale    = estimates.Y_scale;
Y_old      = (Data_nowcast_old.data'-Y_location)./Y_scale;
Y_new      = (Data_nowcast_new.data'-Y_location)./Y_scale;

% Recover dimensions
n   = estimates.dimvec(1);
n_f = estimates.dimvec(2);
p_f = estimates.dimvec(3);
p_e = estimates.dimvec(4);
T   = size(Y_new, 2);

% Find indexes of dates to backcast/nowcast/forecast
t_now    = (find(~isnan(Y_new(i_now, :)), 1, 'last')+3):3:T;
date_GDP = Data_nowcast_new.timekey(t_now);

% Recover parameters and latent variables
param_new        = map_parameter(median(estimates.param_Gibbs, 2), estimates.dimvec);
latent_new       = struct();
latent_new.sigma = mean(estimates.latents.sigmas, 3);
latent_new.s     = mean(estimates.latents.ss, 3);
param_old        = map_parameter(median(estimates_old.param_Gibbs, 2), estimates.dimvec);
latent_old       = struct();
latent_old.sigma = mean(estimates_old.latents.sigmas, 3);
latent_old.s     = mean(estimates_old.latents.ss, 3);

% Update latent variables and construct state-space models
dns_nowcast = NaN(settings.n_GS/settings.n_each, length(t_now));
sigma_old   = NaN(n_f+n, T, settings.n_GS/settings.n_each);
s_old       = NaN(n_f+n, T, settings.n_GS/settings.n_each);
sigma_new   = NaN(n_f+n, T, settings.n_GS/settings.n_each);
s_new       = NaN(n_f+n, T, settings.n_GS/settings.n_each);
del_string  = '';

% Update old latents
rng(321)
for i_draw = 1:settings.n_GS/settings.n_each
    message    = sprintf('Updating volatilities - Draw: %d/%d\n', i_draw, settings.n_GS/settings.n_each);
    fprintf([del_string, message])
    del_string = repmat(sprintf('\b'), 1, length(message));
    latent_old = S_update(param_old, latent_old, Y_old, estimates_old.restrict);
    sigma_old(:, :, i_draw) = latent_old.sigma;
    s_old(:, :, i_draw)     = latent_old.s;
end

% Update new latents
rng(321)
for i_draw = 1:settings.n_GS/settings.n_each
    message    = sprintf('Updating volatilities - Draw: %d/%d\n', i_draw, settings.n_GS/settings.n_each);
    fprintf([del_string, message])
    del_string = repmat(sprintf('\b'), 1, length(message));
    latent_new = S_update(param_new, latent_new, Y_new, estimates.restrict);
    sigma_new(:, :, i_draw) = latent_new.sigma;
    s_new(:, :, i_draw)     = latent_new.s;
end

latent_old.sigma = mean(sigma_old, 3);
latent_old.s     = mean(s_old, 3);
latent_new.sigma = mean(sigma_new, 3);
latent_new.s     = mean(s_new, 3);
SSM_old          = construct_SSM(param_old, latent_old, estimates.restrict);
SSM_new          = construct_SSM(param_new, latent_new, estimates.restrict);

% Compute point nowcast
fprintf('Point nowcast\n')
[nowcast_tmp, forecasts_tmp, news_tmp, weights_tmp] = point_nowcast(Y_old, Y_new, SSM_old, SSM_new, i_now, t_now);
pnt_nowcast   = Y_location(i_now) + Y_scale(i_now)*nowcast_tmp(4, :);
rev_SSM       = Y_scale(i_now)*(nowcast_tmp(2, :) - nowcast_tmp(1, :));
rev_data      = Y_scale(i_now)*(nowcast_tmp(3, :) - nowcast_tmp(2, :));
forecasts     = Y_location + Y_scale.*forecasts_tmp;
news          = Y_scale.*news_tmp;
actual        = news + forecasts;
weights       = Y_scale(i_now)*(weights_tmp./Y_scale);
impacts       = (actual - forecasts) .* weights;

% Compute density nowcast
for i_draw = 1:settings.n_GS/settings.n_each
    message = sprintf('Density nowcast       - Draw: %d/%d\n', i_draw, settings.n_GS/settings.n_each);
    fprintf([del_string, message])
    del_string             = repmat(sprintf('\b'), 1, length(message));
    nowcast_tmp            = density_nowcast(Y_new, SSM_new, i_now, t_now);
    dns_nowcast(i_draw, :) = Y_location(i_now) + Y_scale(i_now)*nowcast_tmp;
end

%% Display nowcast update

fprintf('\n\n\n');
fprintf('Nowcast Update: %s \n', datestr(date_nowcast_new, 'mmmm dd, yyyy'))
fprintf('Nowcast for %s (%s), %s \n\n',spec.SeriesName{i_now},spec.UnitsTransformed{i_now},datestr(date_nowcast_new,'YYYY:QQ'));

% Get result
[T, n]    = size(Data_nowcast_new.data);
releases  = ~isnan(actual(:));
var_idx   = repmat((1:n)', [1, T]);
var_idx   = var_idx(releases);
Forecast = forecasts(releases);
Actual    = actual(releases);
Weight   = weights(releases);
Impact = impacts(releases);

parameter_rev = rev_SSM;
data_rev = rev_data;
revisions_impact = rev_SSM(1) + data_rev(1);
news_table = table(Forecast, Actual, Weight, Impact, 'RowNames', spec.SeriesName(var_idx));


    fprintf('      Impact from parameter and data revisions:      %5.2f\n', revisions_impact)

    fprintf('                     Impact from data releases:      %5.2f\n', sum(news_table.Impact))
    fprintf('                                                  +_________\n')
    fprintf('                                  Total impact:      %5.2f\n\n', ...
            revisions_impact + sum(news_table.Impact))
    fprintf('                  %s nowcast:                  %5.2f\n\n', ...
        datestr(date_nowcast_new, 'mmm dd'), pnt_nowcast(1))
    


fprintf('\n  Impact of Data Releases: \n\n')

disp(news_table)

output = {};
output.news_table = news_table;
output.nowcast = pnt_nowcast(1);
output.date = date_nowcast_new;
save(['output/' fname('Update', date_nowcast_new)], 'output')
