function [h] = plotSD(data, location)

    color = [0 0 0];

    Y = nanstd(data);
    the_mean = nanmean(data);
    % plot a big circle at the mean
    h = plot(location, the_mean, '.', 'markers', 25, 'Color', color);
    plot([location location], [the_mean-Y, the_mean+Y], 'Color', color, 'linewidth', 1);

end