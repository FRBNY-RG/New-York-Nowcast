function visualize_contributions(inputs)

% Extract input for visualization
dimvec            = inputs.dimvec;
i_now             = inputs.i_now;
Data_estimate     = inputs.Data_estimate;
Y_scale           = inputs.Y_scale;
param_Gibbs       = inputs.param_Gibbs;
latents           = inputs.latents;


% Extract dimensions
n   = dimvec(1);
n_f = dimvec(2);

% Define elements for figures
blue   = [ 51,  51, 153]/255;
black  = [  0,   0,   0]/255;
figopt = {'units', 'inches', 'position', [0 0 14 10]};
axopt  = {'box', 'on', 'xgrid', 'on', 'ygrid', 'on', 'fontname', 'palatino', 'fontsize', 14};


% Redefine quantiles for summaries
quant_grid = [0.16, 0.50, 0.84];
n_quant    = length(quant_grid);
mid_quant  = ceil(n_quant/2);

% Extract latent variables
states = latents.states;

% Compute trend, factors, volatilities and outlier indicators
T                 = size(states, 2);
n_each            = n_GS/size(states, 3);
vec_q             = [1, 2, 3, 2, 1]/9;
GDP_factor_draw   = NaN(T, n_GS/n_each);
GDP_factor_q_draw = NaN(T, n_GS/n_each);

for i_draw = 1:n_GS/n_each
    param_draw = map_parameter(param_Gibbs(:, i_draw*n_each), dimvec);
    GDP_factor_draw(:, i_draw) = Y_scale(i_now)*(param_draw.Lambda(i_now, :)*states(1+(1:n_f), :, i_draw));
    for t = 5:T
        GDP_factor_q_draw(t, i_draw) = vec_q*GDP_factor_draw((t-4):t, i_draw);
    end
end
GDP_factor_q = quantile(GDP_factor_q_draw, quant_grid, 2);

%% Factor contribution
% Create figure for factor contribution
close; fig0 = figure();
set(fig0, figopt{:})
ax0   = axes();
fill(ax0, [Data_estimate.timekey; flipud(Data_estimate.timekey)], [GDP_factor_q(:, 1); flipud(GDP_factor_q(:, n_quant))], ...
     blue, 'FaceAlpha', 0.1, 'EdgeColor', 'none') 
hold('on')
plot(ax0, Data_estimate.timekey, GDP_factor_q(:, mid_quant), 'linewidth', 2, 'linestyle', '-', 'color', blue)
hold('on')
miss  = isnan(Data_estimate.data(:, i_now));
plot(ax0, Data_estimate.timekey(~miss), Data_estimate.data(~miss, i_now), 'linewidth', 2, 'linestyle', ':', 'color', black);
hold('off')

% Tune ax handle
set(ax0, axopt{:})
xlim(ax0, [Data_estimate.timekey(1) Data_estimate.timekey(end)])
xlabel(ax0, 'Time')
ylim(ax0, 6*Y_scale(i_now)*[-1, 1])
ylabel(ax0, 'Factor contribution (quarterly)')
legend(ax0, {'', 'Factor contribution', 'GDP growth'}, 'location', 'northwest', 'orientation', 'vertical')

% Save figure
saveas(fig0, [figdir 'decomposition\full_factors_q'], 'png')

% Save figure for recent sample
xlim(ax0, [Data_estimate.timekey(T-240) Data_estimate.timekey(T)])
ylim(ax0, 8*Y_scale(i_now)*[-1, 1])
legend(ax0, {'', 'Factor contribution', 'GDP growth'}, 'location', 'northwest', 'orientation', 'vertical')
saveas(fig0, [figdir 'decomposition\recent_factors_q.png'], 'png')

% Save figure for prepandemic sample
xlim(ax0, [datetime(2000, 1, 1) datetime(2020, 1, 1)])
ylim(ax0, 8*Y_scale(i_now)*[-1, 1])
legend(ax0, {'', 'Factor contribution', 'GDP growth'}, 'location', 'northwest', 'orientation', 'vertical')
saveas(fig0, [figdir 'decomposition\prepandemic_factors_q.png'], 'png')

end
