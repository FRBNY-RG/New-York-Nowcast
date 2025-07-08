function [param, latent] = Gibbs_update(param, latent, Y, prior, restrict)
% GIBBS_UPDATE  Perform Gibbs sampler update to parameters of nowcast
%   model. Model is described in construct_SSM.m and prior distribution
%   in construct_prior.m.
%
%   PARAM = GIBBS_UPDATE(PARAM, LATENT, Y, PRIOR, RESTRICT, UPDATE_PARAM)
%   updates parameters PARAM based on data Y, latent variables LATENT,
%   prior distribution PRIOR and restrictions RESTRICT:
%     PARAM is struct, see construct_SSM.m.
%     LATENT is struct:
%     - LATENT.SIGMA = [SIGMA_F; SIGMA_E] is (N_F+N)xT array with the
%     volatilities of f_t and e_t.
%     - LATENT.S = [S_F; S_E] is (N_F+N)xT array with outlier indicators
%     for f_t and e_t.
%     Y is NxT, contains data.
%     PRIOR is struct, see construct_prior.m.
%     RESTRICT is struct:
%     - RESTRICT.LAMBDA is NxN_F w/NaN if element is unrestricted and
%       the value (usually 0) if restricted.
%     - RESTRICT.PHI is N_FxN_FxP_F w/NaN if element is unrestricted and
%       the value (usually 0) if restricted.
%     - RESTRICT.IOTA is Nx1 with multipliers of the time-varying trend.
%     - RESTRICT.F_ACTIVE is N_FxT w/TRUE in entry (i, t) if factor is
%       active and FALSE otherwise, defaults to TRUE.
%     - RESTRICT.ISQUART is Nx1 w/TRUE if variable is quarterly and false
%       otherwise, defaults to FALSE.
%     UPDATE_PARAM is logical, UPDATE_PARAM=TRUE if parameters are to be
%     updated and UPDATE_PARAM=FALSE otherwise, defaults to TRUE.
%
%   [PARAM, LATENTS] = GIBBS_UPDATE(PARAM, LATENTS, Y, PRIOR, RESTRICT)
%   also returns draw of state variables, volatilities outlier indicators
%   contained in LATENT:
%     LATENT is struct:
%     - LATENT.STATE is N_STATExT.
%     - LATENT.SIGMA = [SIGMA_F; SIGMA_E] is (N_F+N)xT array with the
%     volatilities of f_t and e_t.
%     - LATENT.S = [S_F; S_E] is (N_F+N)xT array with outlier indicators
%     for f_t and e_t.
%
%   Version: 2021 Dec 01 - Matlab R2020a

% Extract parameters from structure
mu      = param.mu;
Lambda  = param.Lambda;
gamma_g = param.gamma_g; %#ok
Phi     = param.Phi;
gamma_f = param.gamma_f;
pi_f    = param.pi_f;
phi     = param.phi;
gamma_e = param.gamma_e;
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

% Define function handle and options for inverse algorithm
symmetrize    = @(A) (A + A')/2;
option        = struct();
option.SYM    = true;
option.POSDEF = true;

% Set support for scale_eps and ps prior
n_s_vals = 100;
s_vals   = [1; linspace(2, 5, n_s_vals-1)'];
s_probs  = [pi_f, (1-pi_f)/(n_s_vals-1)*ones(1, n_s_vals-1); ...
    pi_e, (1-pi_e)/(n_s_vals-1)*ones(1, n_s_vals-1)];

%% LATENT VARIABLES

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
g_t              = state(1, t_skip+(1:T_est));
f_t              = state(n_g_state + (1:n_f), t_skip+(1:T_est));
e_t              = NaN(n, T_est);
e_t(~isquart, :) = state((n_g_state+n_f_state) + (1:(n-n_quart)), t_skip+(1:T_est));
e_t(isquart, :)  = state((n_g_state+n_f_state+n_e_state) + (1:n_quart), t_skip+(1:T_est));

% Extract active-factor indicators
if isfield(restrict, 'f_active')
    F_t = double( restrict.f_active(:, t_skip+(1:T_est)) ) .* f_t;
else
    F_t = f_t;
end

% Extract volatilities and outlier indicators
sigma     = latent.sigma;
s         = latent.s;
sigma_f   = sigma(1:n_f, t_skip+(1:T_est));
s_f       = s(1:n_f, t_skip+(1:T_est));
sigmaXs_f = sigma_f .* s_f;
sigma_e   = sigma(n_f+(1:n), t_skip+(1:T_est));
s_e       = s(n_f+(1:n), t_skip+(1:T_est));
sigmaXs_e = sigma_e .* s_e;

% Construct detrended monthly-equivalent data
y_t = mu + Lambda*F_t + e_t;


%% FACTORS

% Define auxiliary variables
f_lags = NaN(n_f*p_f, T_est-p_f);
for lag = 1:p_f
    f_lags((lag-1)*n_f+(1:n_f), :) = f_t(:, (p_f+1-lag):(T_est-lag));
end
f      = f_t(:, (p_f+1):T_est);

% Update Phi
if (p_f > 0)
    s_Phi        = reshape((1./sigmaXs_f(:, (p_f+1):T_est)), [n_f*(T_est-p_f), 1]);
    r_Phi        = s_Phi .* f(:);
    R_Phi        = s_Phi .* kron(f_lags', eye(n_f));
    vec_Phi      = Phi(:);
    unr_Phi      = isnan(restrict.Phi(:));
    if any(~unr_Phi) % impose restrictions
        r_Phi = r_Phi - R_Phi(:, ~unr_Phi)*vec_Phi(~unr_Phi);
        R_Phi = R_Phi(:, unr_Phi);
    end
    Rr_Phi       = (R_Phi')*r_Phi;
    RR_Phi       = (R_Phi')*R_Phi;
    Pinv_Phi     = symmetrize( linsolve(symmetrize( prior.P_Phi(unr_Phi, unr_Phi) + RR_Phi ), eye(nnz(unr_Phi)), option) );
    vec_m_Phi    = prior.m_Phi(:);
    m_Phi        = Pinv_Phi * (prior.P_Phi(unr_Phi, unr_Phi)*vec_m_Phi(unr_Phi) + Rr_Phi);
    vec_Phi      = mvnrnd(m_Phi, Pinv_Phi)';
    Phi(unr_Phi) = vec_Phi;
    param.Phi    = Phi;
end

% Update sigma_f and s_f
r_f     = f - reshape(Phi, [n_f, n_f*p_f])*f_lags;
for i_f = 1:n_f
    r_f_tmp       = [NaN(1, t_skip+p_f), r_f(i_f, :), NaN(1, T-t_skip-T_est)];
    sigma(i_f, :) = update_vol(r_f_tmp./s(i_f, :), sigma(i_f, :), gamma_f(i_f));
    s(i_f, :)     = update_scl(r_f_tmp./sigma(i_f, :), s_vals, s_probs(i_f, :));
    if i_f == 5
        sigma(i_f, :) = 1; %% TEMP for matfiles_sig1
    end
end
sigma_f = sigma(1:n_f, t_skip+(1:T_est));
s_f     = s(1:n_f, t_skip+(1:T_est));

% Update gamma_f
r_f           = 2*( log(sigma_f(:, (p_f+2):T_est)) - log(sigma_f(:, (p_f+1):(T_est-1))) )';
gamma_f       = update_gam(r_f, prior.nu_f, prior.s2_f);
param.gamma_f = gamma_f;

% Update pi_f
pi_f       = update_ps(s_f', prior.a_f, prior.b_f);
param.pi_f = pi_f;

%% ERRORS

% Define auxiliary variables
e_lags = NaN(n*p_e, T_est-p_e);
for lag = 1:p_e
    e_lags((lag-1)*n+(1:n), :) = e_t(:, (p_e+1-lag):(T_est-lag));
end
e      = e_t(:, (p_e+1):T_est);

% Update phi
if (p_e > 0)
    s_phi     = reshape(1./sigmaXs_e(:, (p_e+1):T_est), [n*(T_est-p_e), 1]);
    r_phi     = s_phi .* e(:);
    R_phi     = zeros(n*(T_est-p_e), n*p_e);
    R_phi(logical(repmat(eye(n), T_est-p_e, p_e))) = e_lags';
    R_phi     = s_phi .* R_phi;
    Rr_phi    = (R_phi')*r_phi;
    RR_phi    = (R_phi')*R_phi;
    Pinv_phi  = symmetrize( linsolve( symmetrize(prior.P_phi + RR_phi), eye(n*p_e), option) );
    m_phi     = Pinv_phi * (prior.P_phi*prior.m_phi(:) + Rr_phi);
    vec_phi   = mvnrnd(m_phi, Pinv_phi)';
    phi(:)    = vec_phi;
    param.phi = phi;
end

% Update sigma_e and s_e
phi_diags = zeros(n, n*p_e);
phi_diags(logical(repmat(eye(n), 1, p_e))) = phi;
r_e       = e - phi_diags*e_lags;
for i = 1:n
    r_e_tmp         = [NaN(1, t_skip+p_e), r_e(i, :), NaN(1, T-t_skip-T_est)];
    sigma(n_f+i, :) = update_vol(r_e_tmp./s(n_f+i, :), sigma(n_f+i, :), gamma_e(i));
    s(n_f+i, :)     = update_scl(r_e_tmp./sigma(n_f+i, :), s_vals, s_probs(n_f+i, :));
end
sigma_e   = sigma(n_f+(1:n), t_skip+(1:T_est));
s_e       = s(n_f+(1:n), t_skip+(1:T_est));
sigmaXs_e = sigma_e .* s_e;

% Update gamma_e
r_e           = 2*( log(sigma_e(:, (p_e+2):T_est)) - log(sigma_e(:, (p_e+1):(T_est-1))) )';
gamma_e       = update_gam(r_e, prior.nu_e, prior.s2_e);
param.gamma_e = gamma_e;

% Update pi_e
pi_e       = update_ps(s_e', prior.a_e, prior.b_e);
param.pi_e = pi_e;


%% MEASUREMENTS

% Define auxiliary variables
y_mu = y_t(:, (p_e+1):T_est) - Lambda*F_t(:, (p_e+1):T_est);
for lag = 1:p_e
    y_mu = y_mu - phi(:, lag) .* ( y_t(:, (p_e+1-lag):(T_est-lag)) - Lambda*F_t(:, (p_e+1-lag):(T_est-lag)) );
end

% Update mu
Rr_mu    = (ones(n, 1)-sum(phi, 2)) .* sum((1./sigmaXs_e(:, (p_e+1):T_est).^2) .* y_mu, 2);
RR_mu    = diag( sum((1./sigmaXs_e(:, (p_e+1):T_est).^2), 2) .* ((ones(n, 1)-sum(phi, 2)).^2) );
Pinv_mu  = linsolve( symmetrize(prior.P_mu + RR_mu), eye(n), option);
m_mu     = Pinv_mu * (prior.P_mu*prior.m_mu + Rr_mu);
mu       = mvnrnd(m_mu, Pinv_mu)';
param.mu = mu;

% Update gamma_g
r_g           = (g_t(:, 2:T_est) - g_t(:, 1:(T_est-1)))';
gamma_g       = update_gam(r_g, prior.nu_g, prior.s2_g);
param.gamma_g = gamma_g;

if (n*n_f < 100) || (n_f == 1)
    % Define auxiliary variables
    y_Lambda = y_t(:, (p_e+1):T_est) - mu;
    f_Lambda = NaN(n*(T_est-p_e), n*n_f, p_e);
    for lag = 1:p_e
        y_Lambda            = y_Lambda - phi(:, lag) .* (y_t(:, (p_e+1-lag):(T_est-lag)) - mu);
        f_Lambda(:, :, lag) = kron(F_t(:, (p_e+1-lag):(T_est-lag))', diag(phi(:, lag)));
    end
    
    % Update Lambda
    s_Lambda           = reshape(1./sigmaXs_e(:, (p_e+1):T_est), [n*(T_est-p_e), 1]);
    r_Lambda           = s_Lambda .* y_Lambda(:);
    R_Lambda           = s_Lambda .* (kron(F_t(:, (p_e+1):T_est)', eye(n)) - sum(f_Lambda, 3));
    vec_Lambda         = Lambda(:);
    unr_Lambda         = isnan(restrict.Lambda(:));
    if any(~unr_Lambda) % impose restrictions
        r_Lambda = r_Lambda - R_Lambda(:, ~unr_Lambda)*vec_Lambda(~unr_Lambda);
        R_Lambda = R_Lambda(:, unr_Lambda);
    end
    Rr_Lambda          = (R_Lambda')*r_Lambda;
    RR_Lambda          = (R_Lambda')*R_Lambda;
    Pinv_Lambda        = symmetrize( linsolve( symmetrize(prior.P_Lambda(unr_Lambda, unr_Lambda) + RR_Lambda), eye(nnz(unr_Lambda)), option) );
    vec_m_Lambda       = prior.m_Lambda(:);
    m_Lambda           = Pinv_Lambda * (prior.P_Lambda(unr_Lambda, unr_Lambda)*vec_m_Lambda(unr_Lambda) + Rr_Lambda);
    vec_Lambda         = mvnrnd(m_Lambda, Pinv_Lambda)';
    Lambda(unr_Lambda) = vec_Lambda;
    param.Lambda       = Lambda;
else
    for i_f = 1:n_f
        % Define auxiliary variables
        i_tmp    = ((1:n_f) ~= i_f);
        y_Lambda = y_t(:, (p_e+1):T_est) - mu - Lambda(:, i_tmp)*F_t(i_tmp, (p_e+1):T_est);
        f_Lambda = NaN(n*(T_est-p_e), n, p_e);
        for lag = 1:p_e
            y_Lambda = y_Lambda - phi(:, lag) .* (y_t(:, (p_e+1-lag):(T_est-lag)) - mu - Lambda(:, i_tmp)*F_t(i_tmp, (p_e+1-lag):(T_est-lag)) );
            f_Lambda(:, :, lag) = kron(F_t(i_f, (p_e+1-lag):(T_est-lag))', diag(phi(:, lag)));
        end
        
        % Update Lambda
        s_Lambda           = reshape(1./sigmaXs_e(:, (p_e+1):T_est), [n*(T_est-p_e), 1]);
        r_Lambda           = s_Lambda .* y_Lambda(:);
        R_Lambda           = s_Lambda .* (kron(F_t(i_f, (p_e+1):T_est)', eye(n)) - sum(f_Lambda, 3));
        vec_Lambda         = Lambda(:, i_f);
        unr_Lambda         = isnan(restrict.Lambda(:, i_f));
        if any(~unr_Lambda) % impose restrictions
            r_Lambda = r_Lambda - R_Lambda(:, ~unr_Lambda)*vec_Lambda(~unr_Lambda);
            R_Lambda = R_Lambda(:, unr_Lambda);
        end
        Rr_Lambda          = (R_Lambda')*r_Lambda;
        RR_Lambda          = (R_Lambda')*R_Lambda;
        P_Lambda           = prior.P_Lambda((i_f-1)*n + (1:n), (i_f-1)*n + (1:n));
        Pinv_Lambda        = symmetrize( linsolve( symmetrize(P_Lambda(unr_Lambda, unr_Lambda) + RR_Lambda), eye(nnz(unr_Lambda)), option) );
        vec_m_Lambda       = prior.m_Lambda(:, i_f);
        m_Lambda           = Pinv_Lambda * (P_Lambda(unr_Lambda, unr_Lambda)*vec_m_Lambda(unr_Lambda) + Rr_Lambda);
        vec_Lambda         = mvnrnd(m_Lambda, Pinv_Lambda)';
        Lambda(unr_Lambda, i_f) = vec_Lambda;
        param.Lambda       = Lambda;
    end
end


%% STORE LATENT VARIABLES

% Clear redundant states
state_clean = [state(1, :); ...
    state(n_g_state+(1:n_f), :); ...
    NaN(n, T)];
state_clean([false(1+n_f, 1); ~isquart], :) = state((n_g_state+n_f_state) + (1:(n-n_quart)), :);
state_clean([false(1+n_f, 1); isquart], :)  = state((n_g_state+n_f_state+n_e_state) + (1:n_quart), :);

% Store latents in struct
latent.state = state_clean;
latent.sigma = sigma;
latent.s     = s;


end