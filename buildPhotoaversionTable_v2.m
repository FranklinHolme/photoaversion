function [summary] = buildPhotoaversionTable_v2(table_name, path, before_stim_length, after_stim_length)    

    % before_stim_length: time to include before the stimulus (in seconds) 
    % after_stim_length: time to after the stimulus (in seconds)
    
    cd(path)
    
    % Summary table with all the mice - used to find spreadsheets of results
    % for each mouse
    opts = detectImportOptions(table_name);
    opts.VariableNamesRange = 'A1';
    summary = readtable(table_name, opts, 'ReadVariableNames', true);
    
    % Manual exclusion of movies
    summary = summary(logical(summary.include), :);
    
    % If there is an autotrack column (indicating whether the automated
    % tracking was successful), exclude the mice for which it wasn't 
    if sum(strcmp('autotrack', summary.Properties.VariableNames)) > 0
       summary = summary(logical(summary.autotrack), :); 
    end

    % Lists the fields of data to consider. Name correspondences between tracking
    % output and summary spreadsheet are implemented by sharing indeces in
    % these lists 
    tracking_fields = {'Theta', 'xT', 'yT', 'X_pixels_', 'Y_pixels_'};
    summary_fields = {'theta', 'xt', 'yt', 'xc', 'yc', 'time', 'xt_filt', 'yt_filt', ...
        'xc_filt', 'yc_filt', 'speedt', 'speedc', 'distance', 'theta_filt', 'theta_filt_wrt_LED', 'theta_filt_wrt_initial', 'theta_cummax'}; 

    % Set up additional columns in the summary table to hold detailed tracking
    % data. They're intialized to the maximum size needed for any data
    % entry, the maximum framerate * the total trial length. 
    for j = 1:length(summary_fields)
        summary.(summary_fields{j}) = NaN .* ones(size(summary, 1), (before_stim_length + after_stim_length) * max(summary.framerate));
    end

    % Iterate through rows of the summary data
    for i = 1:size(summary, 1)

        % Load the tracking results spreadsheet for a trial from this mouse 
        listing_name = [summary.file{i},'.csv'];
        try % In case the tracking results spreadsheet does not exist
            this_mouse = readtable(listing_name);
            disp(listing_name); 
        catch
            error(['Could not read an expected tracking spreadsheet: ', listing_name]);
        end

        % Replace the crazy first entry in the angle data 
        this_mouse.Theta(1) = this_mouse.Theta(2); 

        % Detect movie frames that have low Hough score or where the
        % distance between the tracked position of the hat center and
        % the hole is too long.
        %auto_bad_frames = detectBadFrames(this_mouse);
        auto_bad_frames = []; % disabled for now 

        % For each field of data from the tracking spreadsheets for
        % individual mice
        for j = 1:length(tracking_fields)
            
            % Selects a data field from the tracking spreadsheet
            cur_data = this_mouse.(tracking_fields{j}); 
            
            % Remove the automatically detected bad frames
            cur_data(auto_bad_frames) = NaN;
            
            % Remove manually annotated bad frames 
            try
                if ~isempty(summary.manual_frame_exclusion{i})
                    cur_data(eval(summary.manual_frame_exclusion{i})) = NaN; 
                end
            catch 
                warning('Manual frame exclusion produced an error'); 
            end 
            
            % Moves data from tracking output spreadsheets to the correct row of the summary table.  
            data = selectTrackingData(cur_data, before_stim_length, after_stim_length, summary.light_on_frame(i), summary.framerate(i)); 
            summary.(summary_fields{j})(i,1:length(data')) = data';
            
        end
        
        % Filter the position data. Generates new data columns while
        % keeping the raw, unfiltered data. 
        position_fields = {'xt', 'yt', 'xc', 'yc'};
        filt_position_fields = {'xt_filt', 'yt_filt', 'xc_filt', 'yc_filt'};
        for j = 1:length(position_fields)
            summary.(position_fields{j})(i,:) = summary.(position_fields{j})(i,:) / summary.pixel_to_cm(i); % convert pixels to cm
            summary.(filt_position_fields{j})(i,:) = filterPosition(summary.(position_fields{j})(i,:), summary.framerate(i));
        end 
        
        % Add a time column
        summary.time(i,:) = (1:length(summary.theta(i, :))) / summary.framerate(i); 
        summary.time(i, isnan(summary.theta(i, :))) = NaN;
        
        % Filter the angle data. Generates new columns while keeping the raw, unfiltered data.  
        summary.theta_filt(i,:) = filterAngleData(summary.theta(i, :), summary.framerate(i));
        
        % Zero angle data w.r.t. LED or w.r.t. mouse head angle at light on
        
        % if there is an 'LED_side' column, read it and compare the angle
        % to the appropriate LED reference angle
        if sum(strcmp('LED_side', summary.Properties.VariableNames)) > 0
            if strcmp(summary.LED_side(i), 'L')
                summary.theta_filt_wrt_LED(i,:) = zeroAngleTrace(summary.theta_filt(i,:), 0); % LED on left is 0 degrees 
            else
                summary.theta_filt_wrt_LED(i,:) = zeroAngleTrace(summary.theta_filt(i,:), 180); % LED on right is 180 degrees
            end
        else
            summary.theta_filt_wrt_LED(i,:) = zeroAngleTrace(summary.theta_filt(i,:), 0); % LED is on left when not specified
        end
            
        summary.theta_filt_wrt_initial(i,:) = zeroAngleTrace(summary.theta_filt(i,:), findInitialAngle(summary, i, before_stim_length));
        
        % Calculate angle cumulative maximum: the maximum angle displacement acheived
        % so far at each time point. 
        theta_temp = summary.theta_filt_wrt_initial(i,:);
        theta_temp(1:round(summary.framerate(i)*before_stim_length)) = NaN; % Don't count angle from before light comes on towards distance
        summary.theta_cummax(i,:) = cummax(theta_temp); 
        summary.theta_cummax(i, isnan(summary.theta_filt_wrt_LED(i,:))) = NaN; % 20210215 fix 
        
        % Calculate speed from the hat center and hole positions (in cm/s). 
        summary.speedc(i,:) = positionToSpeed(summary.xc_filt(i,:), summary.yc_filt(i,:), summary.framerate(i));
        summary.speedt(i,:) = positionToSpeed(summary.xt_filt(i,:), summary.yt_filt(i,:), summary.framerate(i));
        
        % Set speed to be the minimum of the speed of the hat center and hole.
        summary.speed(i,:) = min([summary.speedc(i,:); summary.speedt(i,:)]);
        
        % Compute cumulative distance traveled from the speed. dx/dt * dt 
        summary.cumulative_distance(i,:) = NaN * ones(size(summary.speed(i,:))); 
        speed_temp = summary.speed(i,:); 
        speed_temp(1:round(summary.framerate(i)*before_stim_length)) = NaN; % Don't count speed from before light comes on towards distance
        summary.cumulative_distance(i,:) = cumsum(speed_temp ./ summary.framerate(i), 'omitnan');
        
    end
    
    
end

% For each field of data from the raw tracking spreadsheet, takes the
% samples we want to consider and pads missing data with NaNs so that
% all arrays end up being the same length.
function [data_out] = selectTrackingData(data_in, before_stim_length, after_stim_length, light_on_frame, framerate)

    % initialize out array to the correct length
    n_before_out = before_stim_length*framerate;
    n_after_out = after_stim_length*framerate; 
    
    data_out = NaN .* ones(n_before_out + n_after_out, 1);
    
    % Determine what data to select from the input array. Either there is 
    % not enough data (difference between light_on_frame and number of
    % output frames is negative) so start with frame 1 to get data that
    % exists), just enough data (difference is zero, so again start at
    % frame 1), or more than enough (so start at the difference which will
    % be greater than 1). 
    in_start_frame = max(1, light_on_frame - n_before_out); 
    
    % Similar logic to above 
    in_end_frame = min(length(data_in), light_on_frame + n_after_out - 1);
    
    % Select the input data 
    before_data = data_in(in_start_frame:(light_on_frame-1));
    after_data = data_in(light_on_frame:in_end_frame);
    
    % Determine where to put the data in the output array 
    out_start_frame = n_before_out-length(before_data)+1;
    out_end_frame = n_before_out + length(after_data);
    
    % Fill in the data from the input array to the output array
    data_out(out_start_frame:out_end_frame) = [before_data; after_data];     
end

% Zero an angle trace with respect to a reference angle. This takes the
% magnitude of the angular displacement away from the reference while
% discarding the direction
function [theta_zeroed] = zeroAngleTrace(theta, reference) 
    theta_zeroed = NaN * ones(size(theta)); 
    for i = 1:length(theta) 
        theta_zeroed(i) = abs(deltaAngle(theta(i), reference));
    end
end

% Find the initial angle of a mouse's head at the frame before the light
% turned on. Handles the case where the head angle at the exact frame
% before light on is NaN by averaging in time over the immediate vicinity 
function [initial_angle] = findInitialAngle(summary, i, before_stim_length)
    light_on_frame = round(before_stim_length * summary.framerate(i)) + 1;
    initial_angle = summary.theta_filt(i, light_on_frame);
    
    % If the angle at light on is not available, take the median of the first 5
    % frames before the light came on. 
    if isnan(initial_angle)
        initial_angle = nanmedian(summary.theta_filt(i, floor(light_on_frame - summary.framerate(i)/5):light_on_frame)); 
    end
end


