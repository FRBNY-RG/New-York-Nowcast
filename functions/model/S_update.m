function [latent] = S_update(param, latent, Y, restrict)
% S_UPDATE  Update the outlier indicator matrix s
%
%   [LATENT] = S_UPDATE(PARAM, LATENTS, Y, RESTRICT)
%   Returns draw outlier indicators contained in LATENT:
%     LATENT is struct:
%     - LATENT.STATE is N_STATExT; not updated.
%     - LATENT.SIGMA = [SIGMA_F; SIGMA_E] is (N_F+N)xT array with the
%     volatilities of f_t and e_t; not updated.
%     - LATENT.S = [S_F; S_E] is (N_F+N)xT array with updated outlier
%     indicators for f_t and e_t.
%
%   Version: 2023 Mar 22 - Matlab R2020a

% Extract parameters from structure
Lambda  = param.Lambda;
Phi     = param.Phi;
pi_f    = param.pi_f;
phi     = param.phi;
pi_e    = param.pi_e;

% Recover dimensions
T        = size(Y, 2);
[n, n_f] = size(Lambda);
if isempty(Phi), p_f = 0; else, p_f = size(Phi, 3); end
if isempty(phi), p_e = 0; else, p_e = size(phi, 2); end
if isfield(restrict, 'isquart'), isquart = restrict.isquart; else, isquart = false(n, 1); end
n_quart  = nnz(isquart);
t_skip   = p_e + 5*(n_quart>0);
t_est    = ~all(isnan(Y), 1)';
t_est(1:find(t_est, 1, 'last')) = true;
T_est    = nnz(t_est) - t_skip;

% Set support for scale_eps and ps prior
n_s_vals = 100;
s_vals   = [1; linspace(2, 5, n_s_vals-1)'];
s_probs  = [pi_f, (1-pi_f)/(n_s_vals-1)*ones(1, n_s_vals-1); ...
    pi_e, (1-pi_e)/(n_s_vals-1)*ones(1, n_s_vals-1)];

% Form state-space representation and draw states and errors
SSM   = construct_SSM(param, latent, restrict);
state = simulation_smoother(Y, SSM);

% Reconstruct factors and errors
if (n_quart == 0)
    n_g_state = 1;
    n_f_state = max(1, p_f)*n_f;
    n_e_state = max(1, p_e)*n;
else
    n_g_state = 5;
    n_f_state = max(5, p_f)*n_f;
    n_e_state = max(1, p_e)*(n-n_quart);
end
f_t              = state(n_g_state + (1:n_f), t_skip+(1:T_est));
e_t              = NaN(n, T_est);
e_t(~isquart, :) = state((n_g_state+n_f_state) + (1:(n-n_quart)), t_skip+(1:T_est));
e_t(isquart, :)  = state((n_g_state+n_f_state+n_e_state) + (1:n_quart), t_skip+(1:T_est));

% Extract volatilities and outlier indicators
sigma     = latent.sigma;
s         = latent.s;

% Define auxiliary variables
f_lags = NaN(n_f*p_f, T_est-p_f);
for lag = 1:p_f
    f_lags((lag-1)*n_f+(1:n_f), :) = f_t(:, (p_f+1-lag):(T_est-lag));
end
f      = f_t(:, (p_f+1):T_est);

% Update s_f
r_f     = f - reshape(Phi, [n_f, n_f*p_f])*f_lags;
for i_f = 1:n_f
    r_f_tmp       = [NaN(1, t_skip+p_f), r_f(i_f, :), NaN(1, T-t_skip-T_est)];
    s(i_f, :)     = update_scl(r_f_tmp./sigma(i_f, :), s_vals, s_probs(i_f, :));
end

% Define auxiliary variables for errors
e_lags = NaN(n*p_e, T_est-p_e);
for lag = 1:p_e
    e_lags((lag-1)*n+(1:n), :) = e_t(:, (p_e+1-lag):(T_est-lag));
end
e      = e_t(:, (p_e+1):T_est);

% Update s_e
phi_diags = zeros(n, n*p_e);
phi_diags(logical(repmat(eye(n), 1, p_e))) = phi;
r_e       = e - phi_diags*e_lags;
for i = 1:n
    r_e_tmp         = [NaN(1, t_skip+p_e), r_e(i, :), NaN(1, T-t_skip-T_est)];
    s(n_f+i, :)     = update_scl(r_e_tmp./sigma(n_f+i, :), s_vals, s_probs(n_f+i, :));
end

latent.s = s;

end