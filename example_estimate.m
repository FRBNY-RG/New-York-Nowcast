%%% Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script estimates a dynamic factor model (DFM) using a panel of
% monthly and quarterly series.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clear workspace and set paths
close all; clear; clc;
addpath('functions')
addpath(['functions' filesep 'general'])
addpath(['functions' filesep 'model'])
dirs.data = ['data' filesep];

%% User inputs
date_estimate_new = '2023-09-20'; %dataset to use for estimation

date_forecast_old = '2023-12-01'; %pre-allocated space for initvals

%% Load model specification
fname = @(type, d) [type '_' datestr(d, 'yyyy_mm_dd')];

spec = load_spec('model_spec_FRED.csv', 0);

settings = load_settings();

Data_estimate = load([dirs.data fname('Data', date_estimate_new)]);

n = size(spec.SeriesID, 1);

% Extend datasets to cover forecasting period
Data_estimate.timekey    = datetime(Data_estimate.timekey);
date_tmp                 = ((Data_estimate.timekey(end)+calmonths(1)):calmonths:date_forecast_old)'; %#ok
Data_estimate.timekey    = [Data_estimate.timekey; date_tmp];
Data_estimate.data       = [Data_estimate.data; NaN(length(date_tmp), n)];

rng(321)

%% ESTIMATION

fprintf('ESTIMATION\n\n')

initval = load('initval.mat').initval;

% Change location and scale of the data
Y_est      = Data_estimate.data';
Y_location = mean(Y_est(:, Data_estimate.timekey < datetime(2020, 1, 1)), 2, 'omitnan');
Y_scale    = std(Y_est(:, Data_estimate.timekey < datetime(2020, 1, 1)), 0, 2, 'omitnan');
Y_est      = (Y_est-Y_location)./Y_scale;

% Compute dimensions
[n, n_f] = size(spec.Blocks);
p_f      = 4;
p_e      = 1;
dimvec   = [n, n_f, p_f, p_e];
T        = size(Y_est, 2);
t_est    = (Data_estimate.timekey < date_estimate_new);

extend_initval = T -size(initval.sigma, 2);
if extend_initval > 0
    for i = 1:extend_initval
        initval.sigma = [initval.sigma initval.sigma(end)];
        initval.s = [initval.s initval.s(end)];
    end
end 

% Compute indicators for quarterly series
isquart = strcmp(spec.Frequency, 'q');

% Compute restrictions
restrict         = struct();
restrict.Lambda  = spec.Blocks;
restrict.Phi     = NaN(n_f, n_f, p_f);
restrict.iota    = spec.Trend./Y_scale;
restrict.isquart = isquart;

% Prepare input for pandemic factor
i_CoV    = find(strcmpi(spec.BlockNames, 'COVID'));
if (i_CoV > 0)
    t_CoV = or(Data_estimate.timekey < datetime(2020, 3, 1), ...
        Data_estimate.timekey > datetime(2021, 12, 1));
    restrict.Phi(i_CoV, :, :)       = 0;
    restrict.Phi(:, i_CoV, :)       = 0;
    restrict.Phi(i_CoV, i_CoV, :)   = NaN;
    restrict.f_active               = true(n_f, T);
    restrict.f_active(i_CoV, t_CoV) = false;
end

% Set prior
m_Lambda = initval.param.Lambda;
prior    = construct_prior(dimvec, m_Lambda);
prior.P_Phi = prior.P_Phi/5;

% Run Gibbs sampling algorithm
[param_Gibbs, latents] = Gibbs_sampler(Y_est, prior, restrict, initval, settings);


% Store estimates
outputs               = struct();
outputs.spec          = spec;
outputs.dimvec        = dimvec;
outputs.Data_estimate = Data_estimate;
outputs.Y_location    = Y_location;
outputs.Y_scale       = Y_scale;
outputs.isquart       = isquart;
outputs.prior         = prior;
outputs.restrict      = restrict;
outputs.initval       = initval;
outputs.param_Gibbs   = param_Gibbs;
outputs.latents       = latents;
save([dirs.mat fname('Estimates', date_estimate_new)], '-struct', 'outputs')

