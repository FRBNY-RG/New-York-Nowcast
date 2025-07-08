function param = map_parameter(param_vec, dims)
% MAP_PARAMETER  Extract parameter vector of nowcast model into struct. 
% 
%   PARAM = MAP_PARAMETER(PARAM_VEC, DIMS) returns vector assuming param's
%   are vertically concatenated into PARAM_VEC w/dimensions DIMS:
%     PARAM_VEC is N_PARAMx1 array, see vec_parameter.m.
%     DIMS is a vector [N, N_F, P_F, P_E]:
%     - N is number of series (monthly and quarterly).
%     - N_F is number of factors.
%     - P_F is number of lags in VAR model for factors.
%     - P_E is number of lags in AR model for measurement errors.
%     PARAM is struct:
%     - PARAM.mu is Nx1.
%     - PARAM.gamma_g is scalar.
%     - PARAM.Lambda is NxN_F.
%     - PARAM.Phi is N_Fx(N_F*P_F).
%     - PARAM.gamma_f is N_Fx1.
%     - PARAM.pi_f is N_Fx1.
%     - PARAM.phi is NxP_E.
%     - PARAM.gamma_e is Nx1.
%     - PARAM.pi_e is Nx1.
%
%   Version: 2021 Dec 01 - Matlab R2020a

% Recover dimensions
n   = dims(1);
n_f = dims(2);
p_f = dims(3); 
p_e = dims(4);

% Define auxiliary numbers and vectors
n_tmp = cumsum([n, 1, n*n_f, n_f^2*p_f, n_f, n_f, n*p_e, n, n])';
n_vec = [1 + [0; n_tmp(1:(end-1))], n_tmp];

% Map parameters
param         = struct();
param.mu      = reshape(param_vec(n_vec(1, 1):n_vec(1, 2)), [n, 1]);
param.gamma_g = reshape(param_vec(n_vec(2, 1):n_vec(2, 2)), [1, 1]);
param.Lambda  = reshape(param_vec(n_vec(3, 1):n_vec(3, 2)), [n, n_f]); 
param.Phi     = reshape(param_vec(n_vec(4, 1):n_vec(4, 2)), [n_f, n_f, p_f]); 
param.gamma_f = reshape(param_vec(n_vec(5, 1):n_vec(5, 2)), [n_f, 1]);
param.pi_f    = reshape(param_vec(n_vec(6, 1):n_vec(6, 2)), [n_f, 1]);
param.phi     = reshape(param_vec(n_vec(7, 1):n_vec(7, 2)), [n, p_e]);
param.gamma_e = reshape(param_vec(n_vec(8, 1):n_vec(8, 2)), [n, 1]); 
param.pi_e    = reshape(param_vec(n_vec(9, 1):n_vec(9, 2)), [n, 1]); 

end