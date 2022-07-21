function [] = makeFactorPlot(X, h, varargin) 

    % Return nothing and do nothing if the dataset X is empty
    if isempty(X) 
        warning('makeFactorPlot encountered an empty dataset');
        return; 
    end 

    % Select the figure to plot into
    %figure(h); 

    % Data shape: rows (i indeces) are samples, columns (j indeces) are conditions
    numvarargs = length(varargin);
    if numvarargs > 10
        error('makeFactorPlot:TooManyInputs', 'requires at most 6 optional inputs');
    end
    optargs = {false, false, true, true, true, false, 'w', 'k', 5, [0.5, 0.5, 0.5]};
    optargs(1:numvarargs) = varargin; 
    [show_sd, show_mean, show_median, show_lines, show_points, show_boxplot, point_fill_color, point_edge_color, point_size, line_color] = optargs{:};
    
    
    hold on 
    
    if show_lines
        plot(X', 'Color', line_color); 
    end
    
    if show_points
        for j = 1:size(X, 2)

            plot(j* ones(size(X(:,j))), X(:,j), 'o', 'MarkerFaceColor', point_fill_color, 'MarkerEdgeColor', point_edge_color, 'MarkerSize', point_size);

        end
    end 
    
    if show_mean || show_median
        jit = -0.15; 
        for j = 1:size(X, 2)
            if show_mean
                the_stat = nanmean(X(:, j));
            elseif show_median
                the_stat = nanmedian(X(:,j));
            end
           plot([j-0.2, j+0.2], [the_stat, the_stat], 'k', 'LineWidth', 2);
        end
    end
    
    if show_sd %|| show_median
        jit = -0.15; 
        for j = 1:size(X, 2)
           if j == size(X, 2)
               jit = -jit;
           end
           if show_sd
                plotSD(X(:, j), j+jit); 
           elseif show_median
               plotMedian(X(:, j), j+jit);
           end         
        end
    end
            
    
    if show_boxplot
      boxplot(X);
    end
    
    buff = 0.25;
    xticks(1:size(X, 2));
    xlim([1-buff, size(X, 2)+buff]); 
    
    % Formatting
    %set(h, 'Position', [1630, 600, 200+100*(size(X,2)-2), 220]); % add 75 pixels for each additional treatment condition beyond the first two
    set(gca, 'FontSize', 10); 
    box off 
    
end