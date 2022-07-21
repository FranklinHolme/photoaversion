function [h] = plotMedian(data, location)

    color = [0 0 0];

    the_median = nanmedian(data);
    % plot a big circle at the mean
    h = plot(location, the_median, '.', 'markers', 25, 'Color', color);
    
    % Plot lines from the median to the 25th and 75% percentiles 
    %plot([location location], [the_median, the_median], 'Color', color, 'linewidth', 1);

end