function [d, speed_threshold] = calcPhotoaversionTableSummaryStats(d, t_before_stim,...
    t_after_stim, response_window_distance, response_window_turning, n_sds, turn_thresh, instant_angle_time)

    % Set up new columns in the table
    n_movies = size(d, 1);
    d.latency_speed = NaN * ones(n_movies, 1);
    d.latency_movement_bout = NaN * ones(n_movies, 1);
    d.latency_turn = NaN * ones(n_movies, 1); 
    d.latency_distance = NaN * ones(n_movies, 1); 
    d.distance = NaN * ones(n_movies, 1);
    d.max_angle = NaN * ones(n_movies, 1);
    d.manual_latency_observer_1 = NaN * ones(n_movies, 1);
    d.strong_response = NaN * ones(n_movies, 1); 
    d.instant_theta = NaN * ones(n_movies, 1); 
    
    % Calculate speed threshold from dark trials in the dataset
    [dark_sd, dark_median] = getDarkStats(d); 
    speed_threshold = dark_median + dark_sd*n_sds; 
    disp(['Movement speed threshold is ', num2str(speed_threshold), ' cm/s']); 

    for i = 1:n_movies
        
        framerate = d.framerate(i); 
        light_on_frame = round(t_before_stim * framerate) + 1;
        consider_frames_distance = light_on_frame:round(light_on_frame + response_window_distance * framerate - 1); 
        consider_frames_turning = light_on_frame:round(light_on_frame + response_window_turning * framerate - 1); 
        
        % Calculate latency until first fast movement, in seconds.
        % Currently commented out as these take too much time 
        %d.latency_speed(i) = getLatency(d.speed(i, light_on_frame:end), speed_threshold, framerate); % arguments are: trace, threshold, sampling rate
        %d.latency_movement_bout(i) = getMovementBoutLatency(d.speed(i, light_on_frame:end), dark_median, dark_sd, n_sds, framerate);
        
        % Calculate manual latency
        try 
            d.manual_latency_observer_1(i) = (d.observer1_response_frame(i) - d.light_on_frame(i)) / framerate; 
            d.manual_latency_observer_1(i) = min([t_after_stim, d.manual_latency_observer_1(i)], [], 'includenan'); % Floor the latency to max response window 
        catch 
            
            if strcmp('manual_latency_observer_1',d.Properties.VariableNames)
                warning('Calculating manual observer response latency produced an error');
            end 
        end
        
        % Calculate latency until first turn >'turn_thresh' degrees 
        d.latency_turn(i) = getLatency(d.theta_filt_wrt_LED(i, light_on_frame:end), turn_thresh, framerate); 
        
        % Calculate latency until mouse moves more than 1 cm total
        d.latency_distance(i) = getLatency(cumsum(d.speed(i, light_on_frame:end), 'omitnan') / framerate, 1, framerate);
    
        if consider_frames_distance(end) <= length(d.speed(i,:))
            
            % Calculate distance traveled within the response window for distance, in cm 
            d.distance(i) = d.cumulative_distance(i, consider_frames_distance(end));
        end    
            
        if consider_frames_turning(end) <= length(d.theta_cummax(i,:))
            % Calculate maximum angle displacement from initial position,
            % within the response window for turning 
            d.max_angle(i) = nanmax(d.theta_cummax(i,consider_frames_turning(end))); 
        
        else
            
            warning(['There were not enough movie frames to consider a ', num2str(response_window_turning),...
                ' second response window']);
            
        end 
        
        % Decide whether the mouse had a strong response to the light
        % (turned > 90 degrees from LED within 60 seconds) 
        d.strong_response(i) = d.latency_turn(i) <= 60; 
        if isnan(d.strong_response(i)); d.strong_response(i) = false; end
        
        % Measure the instantaneous heading at a particular time 
        d.instant_theta(i) = d.theta_filt_wrt_LED(i, round(instant_angle_time*framerate));
    
    end
end

% Calculate speed threshold using the distribution of speed in all the dark
% trials included in the dataset. Inputs are the data table and the number of 
% standard deviations to set the threshold at 
function [dark_sd, dark_mean] = getDarkStats(d)

    % Find row indices of all the dark trials; 
    is_dark_trial  = d.intensity == 0;
    
    % Extract all the speed data for dark trials
    dark_speed = d.speed(is_dark_trial, :); 
    
    % Reshape the data into a (rather long) 1D array
    dark_speed = dark_speed(:);
    
    % Get the dark stats
    dark_mean = mean(dark_speed, 1, 'omitnan'); 
    dark_sd = std(dark_speed, 1, 'omitnan'); 
    
end


% Calculate latency until a trace exceeds a quantity 
function [latency] = getLatency(trace, threshold, framerate)

    latency = find(trace > threshold, 1, 'first') / framerate; 
    
    if isempty(latency) 
        latency = length(trace) / framerate;
    end

end

% Calculate latency until a movement bout is initiated. First, need to find
% movement bouts, which are when the mouse is moving faster than 3 standard
% deviations of speed in the dark for at least one second. 
function [latency] = getMovementBoutLatency(trace, dark_mean, dark_sd, n_sds, framerate)

    % Define fast movement (conservative) and movement initiation (easy to
    % exceed) thresholds
    fast_move_threshold =  dark_mean + n_sds*dark_sd;
    movement_initiation_threshold = dark_mean; 
    
    % Determine when the mouse is moving fast
    moving_fast = trace > fast_move_threshold;
     
    % Create an array that tracks the mean speed within the following one
    % second. This is to ensure that the movement is prolonged, with a mean
    % speed in the next second above the threshold. 
    movement_bouts = zeros(size(trace));
    for i = 1:length(trace)-framerate
        movement_bouts(i) = mean(trace(i:(i+framerate-1))); 
    end 
    
    % Find the first movement bout by selecting the first index that's part
    % of a prolonged (1 second) period of high speed movement
    first_bout_index = find(movement_bouts >= fast_move_threshold & moving_fast, 1, 'first');
    
    % Find the index where the movement bout started 
    bout_start_index = find(trace(1:first_bout_index) < movement_initiation_threshold, 1, 'last') + 1;
    if isempty(bout_start_index); bout_start_index = first_bout_index; end
    %bout_start_index = first_bout_index;
    
    % Calculate latency
    latency = bout_start_index / framerate;
    
    % If there was never a long enough movement bout, report the full trace
    % length as the latency 
    if isempty(latency) 
        latency = length(trace) / framerate;
    end
    
    if false
        figure
        plot(trace);
        hold on 
        plot(bout_start_index, trace(bout_start_index), '*r');
        trace_temp = trace;
        trace_temp(moving_fast) = NaN;
        plot(trace_temp, '-b');
    end 
        
end
