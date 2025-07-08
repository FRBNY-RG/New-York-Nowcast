function [params, latents] = Gibbs_sampler(Y, prior, restrict, initval, settings)
% GIBBS_SAMPLER  Gibbs sampler for nowcast model. Model is described in
% construct_SSM.m and prior distribution in construct_prior.m.
%
%   PARAMS = GIBBS_SAMPLER(Y, PRIOR, RESTRICT, INITVAL, SETTINGS)
%   returns parameters PARAM based on data Y, prior distribution PRIOR,
%   model restrictions RESTRICT, initial values INITVAL, and numerical
%   settings SETTINGS::
%     Y is NxT, contains data.
%     PRIOR is struct, see construct_prior.m.
%     RESTRICT is struct, see Gibbs_update.m.
%     INITVAL is struct:
%     - INITVAL.PARAM is parameter struct, see map_parameter.m.
%     - INITVAL.SIGMA = [sigma_f; sigma_e] is (N_F+N)xT array with
%       volatilities for f_t and e_t.
%     - INITVAL.S = [s_f; s_e] is (N_F+N)xT array with outlier indicators
%       for f_t and e_t.
%     SETTINGS is struct:
%     - SETTINGS.n_GS is number of updates to be stored.
%     - SETTINGS.N_BURN is number of initial updated burnt-in.
%     - SETTINGS.N_THIN is number of updates between stored.
%     - SETTINGS.N_EACH is number of updates between states/vol's stored.
%     - SETTINGS.PLOT_MCMC is logical, indicates if plot of stored updates
%       is to be displayed, defaults to false.
%     - SETTINGS.plot_each is number of stored updates between plot update.
%     PARAMS is N_PARAMxN_GS array, see vec_parameters.m.
%
%   [PARAMS, LATENTS] = GIBBS_SAMPLER(Y, PRIOR, RESTRICT, INITVAL, SETTINGS)
%   also return draws of latent variables (states, volatilities and outlier
%   indicators):
%     LATENTS is struct:
%     - STATES is N_STATExTx(N_GS/N_EACH).
%     - SIGMAS is (N_F+N)xTx(N_GS/N_EACH).
%     - SS is (N_F+N)xTx(N_GS/N_EACH).
%
%   Version: 2021 Dec 01 - Matlab R2017b

% Determine tasks
need_latents = (nargout > 1);
plot_MCMC    = settings.plot_MCMC;

% Set code parameters and dimensions
n_GS   = settings.n_GS;
n_burn = settings.n_burn;
n_thin = settings.n_thin;
if need_latents, n_each = settings.n_each; end
if plot_MCMC, plot_each = settings.plot_each; end

% Recover dimensions
[n, n_f] = size(initval.param.Lambda);
if isempty(initval.param.Phi), p_f = 0; else, p_f = size(initval.param.Phi, 3); end
if isempty(initval.param.phi), p_e = 0; else, p_e = size(initval.param.phi, 2); end
n_state  = 1 + n_f + n;
T        = size(Y, 2);
n_param  = 1 + n*(1+n_f+p_e+2) + n_f*(n_f*p_f+2);

% Initialize matrices containing parameter draws
params = zeros(n_param, n_GS);
if need_latents
    states = zeros(n_state, T, n_GS/n_each);
    sigmas = zeros(n_f+n, T, n_GS/n_each);
    ss     = zeros(n_f+n, T, n_GS/n_each);
end


%% GIBBS SAMPLING

% Extract initial values
param = initval.param;
sigma = initval.sigma;
s     = initval.s;

% Store latent variables in struct
latent         = struct();
latent.sigma   = sigma;
latent.s       = s;
latent.initsig = sigma;

% Set missing restriction matrices to default values
if ~isfield(restrict, 'Lambda'),   restrict.Lambda   = NaN(n, n_f); end
if ~isfield(restrict, 'Phi'),      restrict.Phi      = NaN(n_f, n_f, p_f); end
if ~isfield(restrict, 'iota'),     restrict.iota     = zeros(n, 1); end
if ~isfield(restrict, 'f_active'), restrict.f_active = true(n_f, T); end

% Initialize MCMC plot
if plot_MCMC
    
    % Define auxiliary indexes for plots
    n_param_M = n*(1+n_f)+1;
    n_param_F = n_f*(n_f*p_f+2);
    n_param_E = n*(p_e+2);
    
    % Set style
    blue             = [ 51, 51, 153]/255;
    red              = [153, 51,  51]/255;
    yellow           = [255, 168, 32]/255;
    green            = [ 51, 153, 51]/255;
    property_M       = cell(n_param_M, 2);
    property_M(:, 1) = {1.5};
    for i = 1:n,                     property_M{i, 2}   = green; end
    property_M{n+1, 2} = yellow;
    for i = (n+1)+(1:n*n_f),         property_M{i, 2}   = red;  end
    property_F       = cell(n_param_F, 2);
    property_F(:, 1) = {1.5};
    for i = 1:(n_f^2*p_f),           property_F{i, 2} = green;  end
    for i = (n_f^2*p_f)+(1:n_f),     property_F{i, 2} = yellow; end
    for i = (n_f^2*p_f+n_f)+(1:n_f), property_F{i, 2} = blue;   end
    property_E       = cell(n_param_E, 2);
    property_E(:, 1) = {1.5};
    for i = 1:(n*p_e),               property_E{i, 2} = green;  end
    for i = (n*p_e)+(1:n),           property_E{i, 2} = yellow; end
    for i = (n*p_e+n)+(1:n),         property_E{i, 2} = blue;   end
    
    % Create subplots for measurement, factor and error parameters
    fig0 = figure();
    set(fig0, 'units', 'normalized', 'position', [0 0.25 0.75 0.95])
    ax_M = subplot(2, 2, 1);
    ax_F = subplot(2, 2, 2);
    ax_E = subplot(2, 2, 3);
    if need_latents
        ax_S = subplot(2, 2, 4);
    end
    drawnow
    
end

del_string = '';
avg_time   = NaN;
tic
for i_GS = (-n_burn):n_GS
    
    % Report progress
    message    = sprintf('Gibbs sampler - Draw: %d/%d - Time left: %d mins %02d secs\n', i_GS, n_GS, floor(avg_time/60), round(mod(avg_time, 60)));
    fprintf([del_string, message])
    del_string = repmat(sprintf('\b'), 1, length(message));
    
    % Run Gibbs-sampling updates
    for i_thin = 1:n_thin
        [param, latent] = Gibbs_update(param, latent, Y, prior, restrict);
        latent.state(6,:) = latent.state(6,:).*restrict.f_active(5,:);
    end
    
    % Store draws
    if (i_GS > 0)
        params(:, i_GS) = vec_parameter(param);
        if need_latents && (mod(i_GS, n_each) == 0)
            states(:, :, i_GS/n_each) = latent.state;
            sigmas(:, :, i_GS/n_each) = latent.sigma;
            ss(:, :, i_GS/n_each)     = latent.s;
        end
    end
    
    % Draw MCMC plot
    if (i_GS > 0) && plot_MCMC && (mod(i_GS, plot_each) == 0)
        plot_M = plot(ax_M, params(1:n_param_M, :)');
        plot_F = plot(ax_F, params(n_param_M+(1:n_param_F), :)');
        plot_E = plot(ax_E, params(n_param_M+n_param_F+(1:n_param_E), :)');
        plot(ax_S, squeeze(mean(ss(1:n_f, :, :) .* sigmas(1:n_f, :, :), 2))', 'linewidth', 2);
        ylabel(ax_M, 'Measurements')
        ylabel(ax_F, 'Factors')
        ylabel(ax_E, 'Errors')
        ylabel(ax_S, 'Average factor volatilities')
        set(plot_M, {'linewidth', 'color'}, property_M)
        set(plot_F, {'linewidth', 'color'}, property_F)
        set(plot_E, {'linewidth', 'color'}, property_E)
        drawnow
    end
    
    % Compute execution time
    avg_time = toc*(n_GS-i_GS)/(i_GS+1+n_burn);
end

% Store latent variables in struct
if need_latents
    latents        = struct();
    latents.states = states;
    latents.sigmas = sigmas;
    latents.ss     = ss;
    latents.initsig = latent.initsig;
end

end