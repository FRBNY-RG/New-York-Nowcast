function nowcast = density_nowcast(Y, SSM, i_now, t_now)
% DENSITY_NOWCAST Perform density nowcast calculations.
% 
%   NOWCAST = DENSITY_NOWCAST(Y, SSM, I_NOW, T_NOW) draws nowcast
%   values NOWCAST based on data Y, state space representation SSM, 
%   cross-sectional index I_NOW and time index T_NOW:
%     NOWCAST is 1xlength(T_NOW), each entry represents a period in T_NOW
%     to be nowcasted/forecasted. 
%     Y is NxT, columns are y_t, missing data is NaN.
%     SSM is struct with state-space matrices, see construct_SSM.m.
%     I_NOW is integer indicating the position of series to be nowcasted.
%     T_NOW is array of integers indicating the periods to be nowcasted or
%     forecasted.
%
%   Version: 2021 Dec 01 - Matlab R2020a

% Determine if measurement equation is time-varying
if (size(SSM.D, 2) > 1), D = SSM.D(i_now, t_now);    else, D = repmat(SSM.D(i_now), [1, length(t_now)]);    end
if (size(SSM.H, 3) > 1), H = SSM.H(i_now, :, t_now); else, H = repmat(SSM.H(i_now, :), [1, 1, length(t_now)]); end

% Draw from nowcast distribution
states  = simulation_smoother(Y, SSM);
nowcast = NaN(1, length(t_now));
for i_t = 1:length(t_now), nowcast(:, i_t) = D(:, i_t) + H(:, :, i_t)*states(:, t_now(i_t)); end

end