function SSM = construct_SSM(param, latent, restrict, var_init)
% CONSTRUCT_SSM  Construct state-space form of nowcast model:
%   Y_(mt) = y_t(~isquarterly),
%   Y_(qt) = (1+2L+3L^2+2L^3+L^4)/9 * y_t(isquarterly),
%   (aggregation is monthly AR to quarterly AR)
%   y_t    = mu + iota*g_t + Lambda*f_t + e_t + o_t,
%   f_t    = Phi * [f_(t-1); ...; f_(t-p_f)] + sigma_(ft) .* s_(ft) .* eps_(ft),
%   e_t    = sum(phi .* [e_(t-1), ..., e_(t-p_e)]) + sigma_(et) .* s_(et) .* eps_(et).
%
%   Time-varying trend and volatilities evolve as
%   g_t               = g_(t-1) + gamma_g * ups_(gt),
%   ln(sigma_(ft).^2) = ln(sigma_(f, t-1).^2) + gamma_f .* ups_(ft),
%   ln(sigma_(et).^2) = ln(sigma_(e, t-1).^2) + gamma_e .* ups_(et).
%
%   Outliers are distributed as
%   s_(ft) ~ Discrete(vals_f, probs_f),
%   s_(et) ~ Discrete(vals_e, probs_e).
%
%   Initial condition is distributed as
%   (f_1, f_0, f_(-1), f_(-2), f_(-3)) ~ N(0_(5*n_f), tau_X*eye(5*n_f)).
%
%   Shocks are distributed as
%   (eps_(ft), eps_(et), ups_(gt), ups_(ft), ups_(et)) ~ N(0_(1+2(n_f+n)), I_(1+2(n_f+n))).
%
%   SSM = CONSTRUCT_SSM(PARAM, LATENT, RESTRICT, VAR_INIT) returns the
%   state-space model matrices based on parameters PARAM, latent variables
%   LATENT, model restrictions RESTRICT and variance of initial condition
%   VAR_INIT:
%     SSM is struct:
%     - SSM.D is Nx1.
%     - SSM.H is NxN_STATE.
%     - SSM.Sigma_eps is NxN.
%     - SSM.F is N_STATExN_STATE.
%     - SSM.G is N_STATExN_SHOCK.
%     - SSM.Sigma_eta is N_SHOCKxN_SHOCK.
%     - SSM.mu_1 is N_STATEx1.
%     - SSM.Sigma_1 is N_STATExN_STATE.
%     PARAM is struct, see map_parameter.m.
%     LATENT is struct:
%     - LATENT.SIGMA = [SIGMA_F; SIGMA_E] is (N_F+N)xT array with the
%     volatilities of f_t and e_t.
%     - LATENT.S = [S_F; S_E] is (N_F+N)xT array with outlier indicators
%     for f_t and e_t.
%     RESTRICT is struct, see Gibbs_update.m.
%     VAR_INIT is the variance of initial value for state variables,
%     defaults to (1e-1)*eye(N_STATE)+(1e2)*BLKD where BLKD is a block
%     diagonal matrix with blocks of ones for each state.
%
%   Version: 2021 Dec 01 - Matlab R2020a

% Extract parameters from structure
mu      = param.mu;
gamma_g = param.gamma_g;
Lambda  = param.Lambda;
Phi     = param.Phi;
gamma_f = param.gamma_f; %#ok
pi_f    = param.pi_f; %#ok
phi     = param.phi;
gamma_e = param.gamma_e; %#ok
pi_e    = param.pi_e; %#ok

% Recover dimensions
[n, n_f] = size(Lambda);
if isempty(Phi), Phi = zeros(n_f); end
p_f      = size(Phi, 3);
if isempty(phi), phi = zeros(n, 1); end
p_e      = size(phi, 2);
if isfield(restrict, 'isquart'), isquart = restrict.isquart; else, isquart = false(n, 1); end
n_quart  = nnz(isquart);

% Compute number of states
if (n_quart == 0)
    n_g_state = 1;
    n_f_state = max(1, p_f)*n_f;
    n_e_state = max(1, p_e)*n;
    n_q_state = 0;
    vec_m     = 1;
    vec_q     = 1;
else
    n_g_state = 5;
    n_f_state = max(5, p_f)*n_f;
    n_e_state = max(1, p_e)*(n-n_quart);
    n_q_state = max(5, p_e)*n_quart;
    vec_m     = [1, 0, 0, 0, 0];
    vec_q     = [1, 2, 3, 2, 1]/9;
end
n_state = n_g_state+n_f_state+n_e_state+n_q_state;

% Extract latent variables
sigma = latent.sigma;
s     = latent.s;
T     = size(sigma, 2);

% Extract trend and active-factor indicators
iota = restrict.iota;
if isfield(restrict, 'f_active')
    f_active = double(restrict.f_active);
else
    f_active = ones(n_f, T);
end
if (n_quart == 0)
    vec_f_active = f_active';
else
    f_always     = all(f_active', 1);
    vec_f_active = [f_active', ...
        [repmat(f_always, [1, 1]); f_active(:, 1:(T-1))'], ...
        [repmat(f_always, [2, 1]); f_active(:, 1:(T-2))'], ...
        [repmat(f_always, [3, 1]); f_active(:, 1:(T-3))'], ...
        [repmat(f_always, [4, 1]); f_active(:, 1:(T-4))']];
end

% Extract volatilities and outlier indicators
sigma_f   = sigma(1:n_f, :);
sigma_e   = sigma(n_f+(1:n), :);
s_f       = s(1:n_f, :);
s_e       = s(n_f+(1:n), :);
sigmaXs_f = sigma_f .* s_f;
sigmaXs_e = sigma_e .* s_e;

% Form state-space matrices for measurement equation
%%% D
D = mu;
%%% H
H = NaN(n, n_state, T);
for t = 1:T
    H(~isquart, :, t) = [kron(vec_m, iota(~isquart)), ...
        kron(vec_m, Lambda(~isquart, :)) .* vec_f_active(t, :), zeros(n-n_quart, n_f_state-length(vec_m)*n_f), ...
        eye(n-n_quart), zeros(n-n_quart, n_e_state-(n-n_quart)), ...
        zeros(n-n_quart, n_q_state)];
    H(isquart, :, t) = [kron(vec_q, iota(isquart)), ...
        kron(vec_q, Lambda(isquart, :)) .* vec_f_active(t, :), zeros(n_quart, n_f_state-length(vec_m)*n_f), ...
        zeros(n_quart, n_e_state), ...
        kron(vec_q, eye(n_quart)), zeros(n_quart, n_q_state-length(vec_q)*n_quart)];
end
%%% Sigma_eps
Sigma_eps = (1e-4)*eye(n);

% Form state-space matrices for transition equation
%%% F
F_g         = [1, zeros(1, n_g_state-1)];
F_f         = [reshape(Phi, [n_f, n_f*p_f]), zeros(n_f, n_f_state-p_f*n_f)];
phi_diags_m = zeros(n-n_quart, n_e_state);
phi_diags_m(logical(repmat(eye(n-n_quart), 1, p_e))) = phi(~isquart, :);
phi_diags_q = zeros(n_quart, n_q_state);
phi_diags_q(logical(repmat(eye(n_quart), 1, p_e))) = phi(isquart, :);
F           = blkdiag([F_g; eye(n_g_state-1), zeros(n_g_state-1, 1)], ...
    [F_f; eye(n_f_state-n_f), zeros(n_f_state-n_f, n_f)], ...
    [phi_diags_m; eye(n_e_state-(n-n_quart)), zeros(n_e_state-(n-n_quart), n-n_quart)], ...
    [phi_diags_q; eye(n_q_state-n_quart), zeros(n_q_state-n_quart, n_quart)]);
%%% G
G = blkdiag([1; zeros(n_g_state-1, 1)], ...
    [eye(n_f); zeros(n_f_state-n_f, n_f)], ...
    [eye(n-n_quart); zeros(n_e_state-(n-n_quart), n-n_quart)], ...
    [eye(n_quart); zeros(n_q_state-n_quart, n_quart)]);
%%% Sigma_eta
Sigma_eta = NaN(1+n_f+n, 1+n_f+n, T-1);
for t = 2:T
    Sigma_eta(:, :, t-1) = blkdiag(gamma_g^2, diag(sigmaXs_f(:, t).^2), ...
        diag(sigmaXs_e(~isquart, t).^2), diag(sigmaXs_e(isquart, t).^2));
end

% Form state-space arrays for initial condition
if (nargin < 4)
    state_group    = cell(1+n_f+n, 1);
    state_group{1} = ones(n_g_state);
    for i_f = 1:n_f,       state_group{1+i_f}             = f_active(i_f, 1) * ones(n_f_state/n_f); end
    for i = 1:(n-n_quart), state_group{1+n_f+i}           = ones(n_e_state/(n-n_quart)); end
    for i = 1:n_quart,     state_group{1+n_f+n-n_quart+i} = ones(n_q_state/n_quart); end
    var_init       = 2*eye(n_state) + 10*blkdiag(state_group{:});
end
mu_1    = zeros(n_state, 1);
Sigma_1 = var_init;

% Pass output as struct
SSM           = struct();
SSM.D         = D;
SSM.H         = H;
SSM.Sigma_eps = Sigma_eps;
SSM.F         = F;
SSM.G         = G;
SSM.Sigma_eta = Sigma_eta;
SSM.mu_1      = mu_1;
SSM.Sigma_1   = Sigma_1;

end