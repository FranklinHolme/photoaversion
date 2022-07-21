function [h] = categoricalLegend(handles, category_data)
% Categories and category data must be strings 

    % Find the first instance of categorey data matching each category 
    categories = unique(category_data, 'stable'); 
    first_category_indices = zeros(1, length(categories));
    for i = 1:length(categories)
         found_indices = strfind(category_data, categories{i}); 
         found_indices = find(not(cellfun('isempty',found_indices)));
         first_category_indices(i) = found_indices(1); 
    end 
    
    % Flip the handles so they match the category data (plotted first to
    % last)
    handles = flip(handles);

    % Create the legend 
    h = legend(handles(first_category_indices), category_data{first_category_indices}, 'Interpreter', 'none'); 
    legend boxoff 
end 