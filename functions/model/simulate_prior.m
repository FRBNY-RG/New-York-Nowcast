function param = simulate_prior(prior, restrict)
% SIMULATE_PRIOR  Draw parameters of nowcast model from prior distribution.
%
%   PARAM = SIMULATE_PRIOR(PRIOR, RESTRICT) returns vector of parameters
%   PARAM from prior distribution PRIOR w/restrictions RESTRICT:
%     PRIOR is struct, see construct_prior.m.
%     RESTRICT is struct, see transition_kernel.m.
%     PARAM is struct, see map_parameter.m.
%
%   Version: 2021 Dec 01 - Matlab R2020a

% Recover dimensions and set function handles
[n, n_f]   = size(prior.m_Lambda);
p_f        = size(prior.m_Phi, 3);
p_e        = size(prior.m_phi, 2);
symmetrize = @(A) full(A + A')/2;

% Draw parameters
mu      = mvnrnd(prior.m_mu(:), symmetrize(pinv(prior.P_mu)))';
gamma_g = 1/sqrt(gamrnd(prior.nu_g/2, 2/(prior.nu_g*prior.s2_g)));
Lambda  = reshape(mvnrnd(prior.m_Lambda(:), symmetrize(pinv(prior.P_Lambda))), [n, n_f]);
Phi     = reshape(mvnrnd(prior.m_Phi(:), symmetrize(pinv(prior.P_Phi))), [n_f, n_f, p_f]);
gamma_f = 1./sqrt(gamrnd(prior.nu_f/2, 2./(prior.nu_f*prior.s2_f), n_f, 1));
pi_f    = betarnd(prior.a_f, prior.b_f, n_f, 1);
phi     = reshape(mvnrnd(prior.m_phi(:), symmetrize(pinv(prior.P_phi))), [n, p_e]);
gamma_e = 1./sqrt(gamrnd(prior.nu_e/2, 2./(prior.nu_e*prior.s2_e), n, 1));
pi_e    = betarnd(prior.a_e, prior.b_e, n, 1);

% Impose restrictions
if isfield(restrict, 'Lambda'), Lambda(~isnan(restrict.Lambda(:))) = restrict.Lambda(~isnan(restrict.Lambda(:))); end
if isfield(restrict, 'Phi'),    Phi(~isnan(restrict.Phi(:)))       = restrict.Phi(~isnan(restrict.Phi(:))); end

% Store parameters
param         = struct();
param.mu      = mu;
param.gamma_g = gamma_g;
param.Lambda  = Lambda;
param.Phi     = Phi;
param.gamma_f = gamma_f;
param.pi_f    = pi_f;
param.phi     = phi;
param.gamma_e = gamma_e;
param.pi_e    = pi_e;

end