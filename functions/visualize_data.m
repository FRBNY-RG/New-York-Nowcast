function visualize_data(inputs)

% Extract input for visualization
figdir           = inputs.figdir;
Data_estimate    = inputs.Data_estimate;
Data_nowcast_old = inputs.Data_nowcast_old;
Data_nowcast_new = inputs.Data_nowcast_new;
spec             = inputs.spec;

% Define elements for figures
red    = [153,  51,  51]/255;
blue   = [ 51,  51, 153]/255;
golden = [218, 165,  32]/255;
figopt = {'units', 'inches', 'position', [0 0 14 10]};
axopt  = {'box', 'on', 'xgrid', 'on', 'ygrid', 'on', 'fontname', 'palatino', 'fontsize', 14};

% Plot series
dates  = sort( datetime(unique([Data_estimate.timekey; Data_nowcast_old.timekey; Data_nowcast_new.timekey])) );
T      = length(dates);
n      = size(Data_estimate.data, 2);
series = NaN(T, n, 3);
series(ismember(dates, Data_estimate.timekey), :, 1)    = Data_estimate.data;
series(ismember(dates, Data_nowcast_old.timekey), :, 2) = Data_nowcast_old.data;
series(ismember(dates, Data_nowcast_new.timekey), :, 3) = Data_nowcast_new.data;

for i = 1:n
    % Create figure
    close all 
    fig0  = figure();
    set(fig0, figopt{:})
    ax0   = axes();
    miss  = any(isnan(series(:, i, 1:3)), 3);
    plot0 = plot(ax0, dates(~miss), squeeze(series(~miss, i, 3:(-1):1)));

    % Tune ax handle
    set(ax0, axopt{:})
    xlim(ax0, [dates(1) dates(end)])
    xlabel(ax0, 'Time')
    ylabel(ax0, sprintf('%s (%s)', spec.SeriesName{i}, spec.UnitsTransformed{i}))
    legend(ax0, {'Current nowcast', 'Previous nowcast', 'Used for estimation'}, 'location', 'best', 'orientation', 'vertical')

    % Tune plot
    set(plot0, 'linewidth', 2)
    set(plot0, {'linestyle'}, {'-'; '-.'; '--'})
    set(plot0, {'color'}, {blue; red; golden})

    % Save figure
    saveas(fig0, [figdir 'full_' num2str(i) '-' Data_estimate.legend(i).mnemonic], 'png')

    % Save figure for recent sample    
    xlim(ax0, [dates(T-240) dates(T)])
    legend(ax0, {'Current nowcast', 'Previous nowcast', 'Used for estimation'}, 'location', 'best', 'orientation', 'vertical')
    saveas(fig0, [figdir 'recent_' num2str(i) '-' Data_estimate.legend(i).mnemonic], 'png')    

    % Save figure for prepandemic sample    
    xlim(ax0, [datetime(2000, 1, 1) datetime(2020, 1, 1)])
    legend(ax0, {'Current nowcast', 'Previous nowcast', 'Used for estimation'}, 'location', 'best', 'orientation', 'vertical')
    saveas(fig0, [figdir 'prepandemic_' num2str(i) '-' Data_estimate.legend(i).mnemonic], 'png')        
end

end