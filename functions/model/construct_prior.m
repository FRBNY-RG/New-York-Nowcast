function prior = construct_prior(dims, m_Lambda)
% CONSTRUCT_PRIOR  Construct default prior for parameters of nowcast model:
%   mu          ~ N(m_mu, P_mu),
%   gamma_g     ~ 1/sqrt( Gamma(nu_g/2, 2./(nu_g*s2_g)) ),
%   vec(Lambda) ~ N(m_Lambda,4 P_Lambda),
%   vec(Phi)    ~ N(m_Phi, P_Phi),
%   gamma_f     ~ 1./sqrt( Gamma_(n_f)(nu_f/2, 2./(nu_f*s2_f)) ),
%   pi_f        ~ Beta_(n_f)(a_f, b_f),
%   vec(phi)    ~ N(m_phi, P_phi),
%   gamma_e     ~ 1./sqrt( Gamma_n(nu_e/2, 2./(nu_e*s2_e)) ),
%   pi_e        ~ Beta_(n_e)(a_e, b_e).
% 
%   PRIOR = CONSTRUCT_PRIOR(DIMS, M_LAMBDA) returns PRIOR based on 
%   dimensions DIMS and prior mean for loadings M_LAMBDA:
%     DIMS is a vector [N, N_F, P_F, P_E]:
%     - N is number of series (monthly and quarterly).
%     - N_F is number of factors.
%     - P_F is number of lags in VAR model for factors.
%     - P_E is number of lags in AR model for measurement errors.
%     M_LAMBDA is NxN_F with prior mean for loadings.
%     PRIOR is struct:
%     - PRIOR.m_mu is Nx1 mean of mu.
%     - PRIOR.P_mu is NxN, precision of mu.
%     - PRIOR.nu_g is degrees of freedom of gamma_g^2.
%     - PRIOR.s2_g is location of gamma_g^2.
%     - PRIOR.m_Lambda is (N*N_F)x1 mean of vec(Lambda).
%     - PRIOR.P_Lambda is (N*N_F)x(N*N_F), precision of vec(Lambda).
%     - PRIOR.m_Phi is (N_F^2*P_F)x1 mean of vec(Phi).
%     - PRIOR.P_Phi is (N_F^2*P_F)x(N_F^2*P_F), precision of vec(Phi).
%     - PRIOR.nu_f is degrees of freedom of gamma_f.^2.
%     - PRIOR.s2_f is N_Fx1 location of gamma_f.^2.
%     - PRIOR.a_f and PRIOR.b_f are prior parameters for pi_f.
%     - PRIOR.m_phi is (N_E*P_E)x1 mean of vec(phi).
%     - PRIOR.P_phi is (N_E*P_E)x(N_E*P_E), precision of vec(phi).
%     - PRIOR.nu_e is degrees of freedom of gamma_e.^2.
%     - PRIOR.s2_e is Nx1 location of gamma_e.^2.
%     - PRIOR.a_e and PRIOR.b_e are prior parameters for pi_e.
% 
%   Version: 2021 Dec 01 - Matlab R2020a

% Recover dimensions
n   = dims(1);
n_f = dims(2);
p_f = dims(3); 
p_e = dims(4);

% Create prior struct
prior = struct();

% Set prior for unconditional means
prior.m_mu = zeros(n, 1);
prior.P_mu = 100*eye(n);

% Set prior for time-varying trend variance
prior.nu_g = 18;
prior.s2_g = 0.0001;

% Set prior for factor loadings
prior.m_Lambda = m_Lambda;
%prior.P_Lambda = 5*eye(n*n_f);
prior.P_Lambda = (10*eye(n*n_f) - (1e-1/n)*kron(eye(n_f), ones(n)));

% Factor VAR coefficients
prior.m_Phi = zeros(n_f, n_f, p_f);
if (p_f > 0), prior.m_Phi(:, :, 1) = eye(n_f); end
Xd_Phi      = [5*diag(kron((1:p_f), ones(1, n_f))); ...
               2*repmat(diag(ones(1, n_f)), [1, p_f]); ...
               2*repmat(ones(1, n_f), [1, p_f])];
P_Phi       = kron((Xd_Phi')*Xd_Phi, eye(n_f));
prior.P_Phi = P_Phi;

% Set prior for factor time-varying volatility variances
prior.nu_f = 2; %6 %2;
prior.s2_f = 0.001; %0.0001 %0.001;

% Set prior for factor outlier probabilities
nper      = 12;
pi_mean   = 1-1/(2*nper);
pi_nobs   = 20; %10*nper; 
prior.a_f = pi_mean*pi_nobs;
prior.b_f = (1-pi_mean)*pi_nobs;

% Set prior for measurement error AR coefficients
prior.m_phi = zeros(n, p_e);
P_phi       = cell(p_e, 1);
for lag = 1:p_e, P_phi{lag} = 25*(lag^2)*eye(n); end
prior.P_phi = blkdiag(P_phi{:});

% Set prior for measurement error time-varying volatility variances
prior.nu_e = 18; % 6 
prior.s2_e = 0.0001;

% Set prior for measurement error outlier probabilities
prior.a_e = pi_mean*pi_nobs;
prior.b_e = (1-pi_mean)*pi_nobs;

% Diffuseness of initial condition
%prior.tau_X = 1000;

end