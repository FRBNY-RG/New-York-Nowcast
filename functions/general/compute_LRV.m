function LRV = compute_LRV(Y, n_lag)
% COMPUTE_LRV Compute LRV using parametric VAR model
%
%   LRV = COMPUTE_LRV(Y, N_LAG) computes long-run variance of data Y using
%   a VAR model with N_LAG lags to approximate the spectrum:
%     Y is NxT, columns are y_t.
%     N_LAG is integer, defaults 0.75*T^(1/3).
%
%   Version: 2020 Dec 28 - Matlab R2017b

[N, T] = size(Y);
if (nargin < 2), n_lag = ceil(0.75*(T^(1/3))); end
y      = Y(:, (n_lag+1):T);
x      = [ones(1, T-n_lag); NaN(N*n_lag, T-n_lag)];
for i_lag = 1:n_lag
    x((2+N*(i_lag-1)):(1+N*i_lag), :) = Y(:, (n_lag+1-i_lag):(T-i_lag));
end
B      = y/x;
U      = y - B*x;
S      = U*(U')/((T-n_lag)-(1+N*n_lag));
IB     = pinv(eye(N) - sum(reshape(B(:, 2:(1+N*n_lag)), [N, N, n_lag]), 3));
LRV    = IB*S*(IB');

end