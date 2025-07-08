function [nowcast, forecasts, news, weights] = point_nowcast(Y_old, Y_new, SSM_old, SSM_new, i_now, t_now)
% POINT_NOWCAST Perform point nowcast calculations.
%
%   NOWCAST = POINT_NOWCAST(Y_OLD, Y_NEW, SSM_OLD, SSM_NEW, I_NOW, T_NOW)
%   updates the nowcast values NOWCAST based on old data Y_OLD, new data
%   Y_NEW, old state space SSM_OLD, new state space SSM_NEW, cross-section
%   index I_NOW and time indexes T_NOW:
%     NOWCAST is 4xlength(T_NOW), each column represents a period in T_NOW
%     to be nowcasted/forecasted. Row 1 has nowcasts/forecasts based on
%     old data with SSM_old, row 2 has nowcasts/forecasts based on old data
%     with SSM_new that isolates the effect of parameter revisions. Row 3
%     has nowcasts/forecasts based on revisions to old data, and row 4 has
%     nowcasts/forecasts based on new data releases.
%     Y_OLD is NxT, columns are y_t, missing data is NaN.
%     Y_NEW is NxT, columns are y_t, missing data is NaN.
%     SSM_OLD and SSM_NEW are struct with state-space matrices, see
%     construct_SSM.m.
%     I_NOW is integer indicating the position of series to be nowcasted.
%     T_NOW is array of integers indicating the periods to be nowcasted or
%     forecasted.
%
%   [NOWCAST, FORECASTS] = POINT_NOWCAST(Y_OLD, Y_NEW, SSM_OLD, SSM_NEW, I_NOW, T_NOW)
%   also computes forecasted values for newly released series FORECASTS:
%     FORECASTS is NxT, each column is forecast of y_t based on old data
%     and revisions, NaN if not newly released.
%
%   [NOWCAST, FORECASTS, NEWS] = POINT_NOWCAST(Y_OLD, Y_NEW, SSM_OLD, SSM_NEW, I_NOW, T_NOW)
%   also computes innovations for newly released series NEWS:
%     NEWS is NxT, each column is innovation to y_t based on old data
%     and revisions, NaN if not newly released.
%
%   [NOWCAST, FORECASTS, NEWS, WEIGHTS] = POINT_NOWCAST(Y_OLD, Y_NEW, SSM, I_NOW, T_NOW)
%   also computes impacts of newly released series on nowcasts/forecasts
%   WEIGHTS:
%     WEIGHTS is NxTxlength(T_NOW), third dimension corresponding to each
%     period in T_NOW.
%
%   Version: 2021 Dec 01 - Matlab R2020a

% Extract dimensions and pre-allocate forecasts, news and weights
[n, T]    = size(Y_old);
forecasts = NaN(n, T);
news      = NaN(n, T);
weights   = NaN(n, T, length(t_now));
nowcast   = NaN(3, length(t_now));

% Determine if measurement equation is time-varying
if (size(SSM_old.D, 2) > 1), D_old = SSM_old.D; else, D_old = repmat(SSM_old.D, [1, T]); end
if (size(SSM_old.H, 3) > 1), H_old = SSM_old.H; else, H_old = repmat(SSM_old.H, [1, 1, T]); end
if (size(SSM_new.D, 2) > 1), D_new = SSM_new.D; else, D_new = repmat(SSM_new.D, [1, T]); end
if (size(SSM_new.H, 3) > 1), H_new = SSM_new.H; else, H_new = repmat(SSM_new.H, [1, 1, T]); end

% Separate revisions from new releases
Y_new               = Y_new(:, 1:length(Y_old));
Y_rev               = Y_new;
Y_rev(isnan(Y_old)) = NaN;
releases            = (~isnan(Y_new) & isnan(Y_old));

% Compute forecasts using old data and old SSM
[~, states]   = fast_smoother(Y_old, SSM_old);
forecasts_tmp = NaN(n, T);
for t = 1:T, forecasts_tmp(:, t) = D_old(:, t) + H_old(:, :, t)*states(:, t); end
nowcast(1, :) = forecasts_tmp(i_now, t_now);

% Compute forecasts using old data and new SSM
[~, states]   = fast_smoother(Y_old, SSM_new);
forecasts_tmp = NaN(n, T);
for t = 1:T, forecasts_tmp(:, t) = D_new(:, t) + H_new(:, :, t)*states(:, t); end
nowcast(2, :) = forecasts_tmp(i_now, t_now);

% Compute forecasts using revised data and new SSM
[~, states]   = fast_smoother(Y_rev, SSM_new);
forecasts_tmp = NaN(n, T);
for t = 1:T, forecasts_tmp(:, t) = D_new(:, t) + H_new(:, :, t)*states(:, t); end
nowcast(3, :) = forecasts_tmp(i_now, t_now);

% Store forecasts and news
forecasts(releases) = forecasts_tmp(releases);
news(releases)      = Y_new(releases) - forecasts(releases);

% Compute forecasts using new data
[~, states]   = fast_smoother(Y_new, SSM_new);
forecasts_tmp = NaN(n, T);
for t = 1:T, forecasts_tmp(:, t) = D_new(:, t) + H_new(:, :, t)*states(:, t); end
nowcast(4, :) = forecasts_tmp(i_now, t_now);

% Compute weights
if (nargout > 3)
    for i = 1:n
        for t = 1:T
            if releases(i, t)
                
                % Construct dummy matrix for difference smoothing
                Y_dummy               = zeros(size(Y_new));
                Y_dummy(i, t)         = 1;
                Y_dummy(isnan(Y_new)) = NaN;
                Y_dummy               = Y_dummy + D_new;
                
                % Run smoother to recover weights
                [~, states]      = fast_smoother(Y_dummy, SSM_new);
                for i_t = 1:length(t_now)
                    weights(i, t, i_t) = H_new(i_now, :, t_now(i_t))*states(:, t_now(i_t));
                end
                
            end
        end
    end
end

end