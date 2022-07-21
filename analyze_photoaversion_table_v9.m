t_before_stim = 30; % Time to include before light stimulus, in seconds 
t_after_stim = 90; % Time to include after light stimulus, in seconds

% Plots to make
polar_plots = true;

% Moment at which to measure all the statistics 
measure_moment = 15; % s 

% Plotting limits (s) 
distance_limit = 90; 
angle_limit = 15; 

% Set the x limit on distance and angle trace plots 
xlimit = [0 90]; % sec. 

n_sds = 2; % Number of standard deviations past the mean for a speed to count as a 'response' movement 
turn_thresh = 90; % Threshold for a strong response 

% A dedicated intensity dependence dataset? (one trial per mouse, except for dark) 
intensity_dependence = false; 

% Compare genotypes? (each mouse tested at multiple light intensities)
plot_genotypes = true; 

% Opn4 Cre/+ Brn3b DTA +/wt or wt/wt
opn4_cre_brn3b_dta = true;

% Examine TRPC 3 6 7? 
TRPC3_6_7 = false; 

% Examine melanopsin ko? 
melanopsin_ko = false;          

% Enucleated? 
enuc = false;

% Examine connexin ko? 
connexin_45 = false; 
connexin_30point2 = false;

% Or compare eye injections? 
injections = false; 

% Select an eye injection? 
injection_select = false;
injection = 'ctrl';

% Select an experiment type?
select_experiment_type = false;
experiment_type = 'ctrl_mfa';

% Compare some other criterion (e.g. chamber types)
arbitrary_selection = false; 



%% Apply primary selections to the data 
cd('D:\');
path = 'general code\photoaversion';

if opn4_cre_brn3b_dta
    data_guide_name = 'Opn4Cre Brn3b DTA.xlsx';
    genotype_1 = 'Opn4_Cre/+; Brn3b_+/+';
    genotype_2 = 'Opn4_Cre/+; Brn3b_DTA/+';
    genotype_3 = '';
end

% Intensity dependence 
if intensity_dependence
    data_guide_name = 'intensity dependence.xlsx'; 
end 

% Enucleated
if enuc
    data_guide_name = 'enuc.xlsx';
end

% TRPC 3 6 7 
if TRPC3_6_7
    data_guide_name = 'TRPC_3_6_7.xlsx';
end 

% Eye injection path
if injections 
    data_guide_name = 'SCH_MFA.xlsx';
end 

% Melanopsin ko path 
if melanopsin_ko
    genotype_1 = 'Opn4 Cre/+';
    genotype_2 = 'Opn4 Cre/Cre';
    genotype_3 = 'none';
    data_guide_name = 'Opn4.xlsx';
end

% Connexin ko path 
if connexin_45 
    genotype_1 = 'Opn4 Cre/+; Cx45 fl/+';
    genotype_2 = 'Opn4 Cre/+; Cx45 fl/fl';
    genotype_3 = 'Opn4 Cre/+; Cx45 +/+';
    data_guide_name = 'Cx45.xlsx';
end

if connexin_30point2
    genotype_1 = 'wt';
    genotype_2 = 'ko';
    genotype_3 = 'none';
    data_guide_name = 'Cx30point2.xlsx';
end

% Make the photoaversion table from the summary data sheet and the
% tracking spreadsheets from each mouse 
d_orig = buildPhotoaversionTable_v2(data_guide_name, path, t_before_stim, t_after_stim);

% Date
%d_orig = d_orig(d_orig.date == 20211109, :);

% Exclude entries where the intensity has not been filled in 
d_orig = d_orig(~isnan(d_orig.intensity), :); 

% Convert irradiance from /cm^2 to /um^2
%d_orig.intensity = d_orig.intensity - 8;
%d_orig.intensity(d_orig.intensity < 0) = 0;

% Include only entries from after the LED was collimated
d_orig = d_orig(strcmp(d_orig.LED_optics, 'collimated'), :);

% Calculate a set of summary statistics for each row of the table and add
% the statistics to new rows. 
[d_orig, movement_speed_threshold] = calcPhotoaversionTableSummaryStats(d_orig, ...
    t_before_stim, t_after_stim, measure_moment, measure_moment, n_sds, turn_thresh, measure_moment);

if intensity_dependence
    d = d_orig(d_orig.trial == 1 | d_orig.intensity == 0, :); % only include single trials of light stimuli 
    d = d(mod(d.intensity, 1) == 0, :); % Only Log intensities that are integers 
    d = d(d.date >= 200609, :);
    for i = 1:size(d, 1); d.distribution{i} = 'same'; end
end 

% Round intensities to integer values
d_orig.intensity = round(d_orig.intensity);

% Select genotypes 
if plot_genotypes
    %d = d_orig(strcmp(d_orig.genotype, genotype_1) | strcmp(d_orig.genotype, genotype_2), :);
    d = d_orig; %bypass this selection step (200210) 
    %d = d_orig(d_orig.date == 2021109, :); % Temporary (20200402) 
end 

if enuc 
    d = d_orig;
end 

% Select data in some arbitrary way
if arbitrary_selection 
    d = d_orig(d_orig.date == 20211107, :);
end 

% Select injection
if injection_select
    d = d_orig(strcmp(d_orig.injection, injection), :);
elseif injections
    d = d_orig(d_orig.intensity ~= 13, :); % Don't need the 13 log unit stim in SCH experiments 
end 



%% Apply secondary selections to the data

% Select a type of experiment 
if select_experiment_type
    d = d(strcmp(d.experiment_type, experiment_type), :);
end 

% For Cx30.2 experiments, illumination of the behavior chamber wasn't set
% up well enough after 170623. Also, het's aren't needed.
if connexin_30point2
    d = d(d.date <= 170623, :);
    %d = d(d.date == 170623, :);
    d = d(~strcmp(d.genotype, 'het'), :);
end 

% For Cx45 experiments, we need one copy of Cre to control for copy number.
% Don't need 13 log unit stims as these are below threshold
if connexin_45
    d = d(~strcmp(d.genotype, 'Opn4 +/+; Cx45 fl/fl'), :);
    d = d(d.intensity ~= 13, :);
end

%% Set up plotting color codes for the data categories 

% Set up default color code
for i = 1:length(d.genotype)
    d.color_code{i} = 'k';
end 
    
% Set up color code for injections
if injections
    for i = 1:length(d.injection)
        cur_manip = d.injection{i};
        switch cur_manip
            case 'ctrl'
                inj_color = 'k';
            case 'mfa'
                inj_color = 'r';
            case 'sch'
                inj_color = 'c';
            otherwise
                inj_color = 'g';
                warning(['Unexpected injection descriptor: ', cur_manip]);
        end
        d.color_code{i} = inj_color;
        d.distribution{i} = cur_manip;
    end
end

% Set up color code for manipulations 
if enuc 
    for i = 1:length(d.manipulation)
        cur_manip = d.manipulation{i};
        switch cur_manip
            case 'sham'
                manip_color = 'k';
            case 'bilat'
                manip_color = 'r';
            case 'unilat'
                manip_color = 'b';
            otherwise
                manip_color = 'g';
                warning(['Unexpected manipulation descriptor: ', cur_manip]);
        end
        d.color_code{i} = manip_color;
        d.distribution{i} = cur_manip; 
    end
end

%% Count missing alleles from the genotype metadata. Set up the color codes for plots
if plot_genotypes && ~injection_select
    
    for i = 1:length(d.genotype)
        
        % Look at allele columns of table
        if TRPC3_6_7
            
            % Color code by number of missing TrpC6 and TrpC7 alleles.
            % Disregard TrpC3 for now but remain aware it could confer
            % motion coordination deficit and check this. 
            n_missing_alleles = 0; 
            
            switch d.trpc6{i}
                case 'wt'
                    % add nothing
                case 'het'
                    n_missing_alleles = n_missing_alleles + 1; 
                case 'hom'
                    n_missing_alleles = n_missing_alleles + 2;
            end 
            
            switch d.trpc7{i}
                case 'wt'
                    % add nothing
                case 'het'
                    n_missing_alleles = n_missing_alleles + 1; 
                case 'hom'
                    n_missing_alleles = n_missing_alleles + 2;
            end 
            
            d.n_missing_alleles{i} = num2str(n_missing_alleles); 
            
            switch n_missing_alleles
                case 0
                    genotype_color = 'k';
                case 1 
                    genotype_color = 'b';
                case 2
                    genotype_color = 'g';
                case 3 
                    genotype_color = [240, 90, 40] ./ 255;
                case 4 
                    genotype_color = 'r';
            end 
            
            d.distribution{i} = num2str(n_missing_alleles); 
                 
            
        else
            switch d.genotype{i}
                case genotype_1
                    genotype_color = 'k';
                case genotype_2
                    genotype_color = 'r';
                case genotype_3 
                    genotype_color = 'b';
                otherwise
                    genotype_color = 'g';
                    warning(['Unexpected genotype: ', d.genotype{i}]);
                     
            end
            
            % if experimental distribution not specified, assume it's genotype
            if sum(strcmp('distribution', d.Properties.VariableNames)) == 0 
                d.distribution{i} = d.genotype{i};
            end 
            
        end
        
        d.color_code{i} = genotype_color;
    end
end

% For TrpC 3 6 7 mice, we only want to compare TrpC6/7 KOs to WTs. Remove
% this if you want to examine hets. Also sort the table by genotype so that
% the KOs plot on top 
if TRPC3_6_7
    d = d(strcmp(d.n_missing_alleles, '0') | strcmp(d.n_missing_alleles, '4'), :); 
    [~, sort_idx] = sort(d.genotype);
    reverse_sort_idx = abs(sort_idx - (max(sort_idx) + 1)); 
    d = d(reverse_sort_idx, :); 
end

%% Reject all trials from mice that turned >90 degrees from LED in <= 60 s in the dark
% or were more than 45 degrees from LED at trial onset 
% Except for enucleated mice - there the priority is minimizing # of
% animals needed 
illegal_mice = [];
for i = 1:length(d.genotype)
    if d.intensity(i) == 0 && ~intensity_dependence && ~enuc && opn4_cre_brn3b_dta% because intensity dependence dataset has only one trial per mouse
        if sum(d.theta_filt_wrt_LED(i, d.light_on_frame(i):round(d.framerate(i) * 60 + d.light_on_frame(i))) > 90) > 0 ... 
            || d.theta_filt_wrt_LED(i, d.light_on_frame(i)) > 45
            illegal_mice = [illegal_mice, d.mouse_id(i)];
        end
    end
end 

for i = 1:length(illegal_mice) % I wanted to do this all in a while loop, but somehow it didn't work. 
    d = d(d.mouse_id ~= illegal_mice(i), :); 
    disp(['Rejected mouse ', num2str(illegal_mice(i)), ' for turning >90 degrees in the first 60 s.']);
end

% Set up color code for arbitrary selection 
if arbitrary_selection
    for i = 1:length(d.chamber)
        switch d.chamber{i}
            case 'all_clear'
                chamber_color = 'r';
            case 'front_clear'
                chamber_color = 'k';
            otherwise
                chamber_color = 'g';
                warning('Unexpected chamber type'); 
        end 
        d.color_code{i} = chamber_color;
    end 
end 

%% Collect some basic facts
mouse_ids = unique(d.mouse_id, 'stable');
intensities = sort(unique(d.intensity, 'stable'));
n_intensities = length(intensities);
n_ids = length(mouse_ids);
unique_distributions = unique(d.distribution);
%unique_colors = unique(d.color_code);

%% Plots to 1) illustrate the effect of filtering on head angle data and 2) display speed and angle displacement
if polar_plots     
    for i = 1:size(d, 1)

        % Generate a polar plot of raw and filtered head angle
        hpolar = figure;
        head_angle_radians = d.theta(i, :) * pi/180;
        head_angle_radians_filt = d.theta_filt(i, :) * pi/180;
        time = (1:length(head_angle_radians)) / d.framerate(i);
        polarplot(head_angle_radians, time, 'r');
        hold on 
        polarplot(head_angle_radians_filt, time, 'k', 'LineWidth', 1);
        set(gcf, 'Position', [0        401         372         190]);
        disp(d.file{i});
        rticks(30:30:90);
        
        rticklabels({'0 sec.', '30 sec.'}); 
        thetaticks(0:90:270);
        thetaticklabels({['0', char(176)], ['90', char(176)], ['180', char(176)], ['270', char(176)]});
        set(gca, 'FontSize', 8);
        
        if d.intensity(i) == 0
            title('dark');
        else
            title([num2str(d.intensity(i)), ' log units']);
        end
       
        hold off 
%         hold on
%         plot(d.xc(i,:), d.yc(i,:), '.-r');
%         plot(d.xc_filt(i,:), d.yc_filt(i,:), '.-k');
%         xlim([0 12.5]); % Chamber dimensions are 12.5 x 3.5 cm
%         ylim([0 3.5]); 
%         aspect_ratio = 12.5/3.5;
%         set(gcf, 'Position', [115  175   200*aspect_ratio   200]);
%         hold off
%         xlabel('X head position (cm)');
%         ylabel('Y head position (cm)'); 
        
        % Plot head speed 
%         hspeed = figure;
%         plot(d.speed(i,:), '.-k');
%         set(gcf, 'Position', [885  810   1035   175]);
%         ylim([0 10]);
%         formatPhotoaversionTrace(t_before_stim, d.framerate(i))
%         ylabel('Head speed (cm/s)');
        
        % Plot cumulative distance traveled 
%         hdistance = figure;
%         distance = cumsum(d.speed(i,:), 'omitnan') / d.framerate(i); 
%         plot(distance, '.-k');
%         set(gcf, 'Position', [885   180   1020   215]);
%         formatPhotoaversionTrace(t_before_stim, d.framerate(i))
%         ylabel('Distance traveled (cm)'); 
        
        % Plot head angle displacement
        hangle = figure;
        plot(time - t_before_stim, d.theta_filt_wrt_initial(i,:), '-k');
        hold on 
        %plot(time - t_before_stim, d.theta_filt_wrt_LED(i,:), '-b');
        plot(time - t_before_stim, d.theta_cummax(i,:), '--k');
        set(gcf, 'Position', [371   400   500   190]);
        ylim([-5 180]); 
        yticks([0 90 180]);
        xlim([-30 90]);
        xticks(-30:30:90);
        %formatPhotoaversionTrace(t_before_stim, d.framerate(i))
        ylabel(['Head angle (', char(176), ')']); 
        xlabel('Time (sec. from light onset)'); 
        %legend({'from t = 0', 'from LED', 'cumulative max.'}, 'Location', 'northwest'); 
        legend({'from t = 0', 'cumulative max.'}, 'Location', 'northwest'); 
        box off 
        set(gca, 'TickDir', 'out', 'FontSize', 8); 
        title(d.file(i), 'Interpreter', 'none');
        
        pause();
        close(hpolar);
        %close(hposition);
        %close(hspeed);
        close(hangle);
        %close(hdistance); 
        
    end
end



%% Plot combined data 

% Only include the first trial (Could cause some problems where trial
% notation is different, e.g. a new trial number for each irradiance)
d = d(d.trial == 1, :);

% Plot cumulative maximum head angle, computed from the moment of light on, for each intensity
fig = figure; 
set(fig,'defaultLegendAutoUpdate','off');
h_angle = tiledlayout(1,n_intensities, 'TileSpacing', 'compact'); 
h_angle.XLabel.String = 'Time (sec. from light onset)';
h_angle.YLabel.String = ['Max. turn (', char(176), ')'];
h_angle.XLabel.FontSize = 10;
h_angle.YLabel.FontSize = 10;

for intensity_idx = 1:n_intensities
    intensity = intensities(intensity_idx);
    nexttile;
    
    ylim([0 180]); 
    h_vline = vline(measure_moment);
    h_vline.Color = 'k';
    hold on 
    
    angle = d.theta_cummax(d.intensity == intensity, :);
    framerate = d.framerate(d.intensity == intensity);
    color_code = d.color_code(d.intensity == intensity);
    genotype = d.genotype(d.intensity == intensity);
    distributions = d.distribution(d.intensity == intensity);
    max_angle = NaN(size(distributions));
    
        
    for i = 1:size(angle, 1)
        cur_angle = angle(i, :);
        if sum(isnan(cur_angle)) ~= length(cur_angle) 
            cur_angle = cur_angle(~isnan(cur_angle));
            time = (1:length(cur_angle)) / framerate(i);
            max_angle(i) = min(cur_angle(time > 15));
            h_line = plot(time, cur_angle, 'LineWidth', 0.5);
            hold on 
            set(h_line, {'color'}, color_code(i)); 
        end
    end
    
    %title(['10^{', num2str(intensity), '} photons/cm^2/s']);
%     if intensity == 0
%         title('Dark');
%     end
    yticks([0, 90, 180]);
    yticklabels({'0', '90', '180'});
    ylim([-5, 180]);
    %h = vline(measure_moment);
    %h.Color = [0.5 0.5 0.5];
    xticks(0:30:90);
    set(gca, 'FontSize', 8, 'TickDir', 'out'); 
    xlim(xlimit); 
    box off 
    h = findobj(gca,'Type','line');
    if intensity_idx == 1
        if plot_genotypes
            if TRPC3_6_7
                h_d_legend = categoricalLegend(h, d.n_missing_alleles(d.intensity == intensity));
                h_d_legend.String = {'genetic ctrl.', 'TrpC 6&7 KO'}; 
                h_d_legend.FontSize = 8;
            else
                categoricalLegend(h, d.genotype(d.intensity == intensity)); 
            end
        end
        if injections 
            categoricalLegend(h, d.injection(d.intensity == intensity)) 
        end
    else % Get rid of y axis 
%         ax = gca; 
%         ax.YAxis.Visible = 'off';
    end
    
    % Move tick labels closer to ticks
    a=gca;
    a.XRuler.TickLabelGapOffset = -1;    % negative numbers move the ticklabels down (positive -> down)
    a.YRuler.TickLabelGapOffset = 1;    % negative numbers move the ticklabels right (negative -> left)

    % Add mean +- 95% CI 
%     for i = 1:length(unique_distributions)
%         this_max_angle = max_angle(strcmp(unique_distributions{i}, distributions));
%         plot(15, nanmean(this_max_angle),...
%             '.', 'MarkerSize', 16, 'Color', unique_colors{i});
%         
%         % Bootstrapped 95% confidence intervals around mean 
%         if ~isempty(this_max_angle)
%             cis = bootci(1000, @mean, this_max_angle);
%             plot([15, 15], cis, 'Color', unique_colors{i}, 'LineWidth', 1);
%         end 
%         
%     end
    hold off 
    
    if intensity == 0
        title('dark');
    else
        title([num2str(intensity), ' log units']);
    end 
    
end
set(gcf, 'Position', [0       642         937         198]); % For three panels, resize width for other # of panels to match
set(gcf, 'Renderer', 'painters');


% Plot head angle from LED, computed from the moment of light on, for each intensity
% Overlay the different distributions for the same intensity
fig = figure; 
set(fig,'defaultLegendAutoUpdate','off');
h_angle = tiledlayout(1,n_intensities, 'TileSpacing', 'compact'); 
h_angle.XLabel.String = 'Time (sec. from light onset)';
h_angle.YLabel.String = ['Head angle (', char(176), ')'];
h_angle.XLabel.FontSize = 10;
h_angle.YLabel.FontSize = 10;

for intensity_idx = 1:n_intensities
    intensity = intensities(intensity_idx);
    nexttile;
    
    angle = d.theta_filt_wrt_LED(d.intensity == intensity, :);
    time = d.time(d.intensity == intensity, :);
    framerate = d.framerate(d.intensity == intensity);
    color_code = d.color_code(d.intensity == intensity);
    genotype = d.genotype(d.intensity == intensity);
    distributions = d.distribution(d.intensity == intensity);
    max_angle = NaN(size(distributions));
        
    for i = 1:size(angle, 1)
        cur_angle = angle(i, :);
        if sum(isnan(cur_angle)) ~= length(cur_angle) 
            %cur_angle = cur_angle(~isnan(cur_angle));
            time = (1:length(cur_angle)) / framerate(i);
            max_angle(i) = min(cur_angle(time > 15));
            h_line = plot(time - t_before_stim, cur_angle, 'LineWidth', 0.5);
            hold on 
            set(h_line, {'color'}, color_code(i)); 
        end 
    end
    
    %title(['10^{', num2str(intensity), '} photons/cm^2/s']);
%     if intensity == 0
%         title('Dark');
%     end
    ylim([0 180]); 
    yticks([0, 90, 180]);
    yticklabels({'0', '90', '180'});
    ylim([-5, 180]);
    %h = vline(measure_moment);
    %h.Color = [0.5 0.5 0.5];
    xticks(0:30:90);
    set(gca, 'FontSize', 8, 'TickDir', 'out'); 
    xlim(xlimit); 
    box off 
    h = findobj(gca,'Type','line');
    if intensity_idx == 1
        if plot_genotypes
            if TRPC3_6_7
                h_d_legend = categoricalLegend(h, d.n_missing_alleles(d.intensity == intensity));
                h_d_legend.String = {'genetic ctrl.', 'TrpC 6&7 KO'}; 
                h_d_legend.FontSize = 8;
            else
                
                
                categoricalLegend(h, d.genotype(d.intensity == intensity)); 
            end
        end
        if injections 
            categoricalLegend(h, d.injection(d.intensity == intensity)) 
        end
    else % Get rid of y axis 
        ax = gca; 
        ax.YAxis.Visible = 'off';
    end
    
    % Move tick labels closer to ticks
    a=gca;
    a.XRuler.TickLabelGapOffset = -1;    % negative numbers move the ticklabels down (positive -> down)
    a.YRuler.TickLabelGapOffset = 1;    % negative numbers move the ticklabels right (negative -> left)

    hold off 
    
    title(intensity);
    
end
set(gcf, 'Position', [0       568         667         190]); % For three panels, resize width for other # of panels to match
set(gcf, 'Renderer', 'painters');



%% Plot raw data 

% Plot mouse position tracks and head angle displacements separately for
% each mouse, at each intensity. 
if false
    h_position = figure;
    h_angle = figure; 
    set(h_angle, 'Renderer', 'painters', 'Position', [1240 558 560 420]);
    set(h_position, 'Renderer', 'painters', 'Position', [680 795 560 183]); 
    
    for id_idx = 1:n_ids
       
        for intensity_idx = 1:n_intensities             
                
            mouse_id = mouse_ids(id_idx);
            intensity = intensities(intensity_idx);
            selector = d.mouse_id == mouse_id & d.intensity == intensity; 

            if sum(selector) > 0
                % Plot mouse position track
                figure(h_position);
                plot(d.xc_filt(selector, :), d.yc_filt(selector, :), d.color_code{selector});
                xlim([0 12.5]); % Chamber dimensions are 12.5 x 3.5 cm
                ylim([0 3.5]); 
                
                % Make time trace
                time = (1:length(d.theta_filt_wrt_LED(selector, :))) / d.framerate(selector); 

                % Plot mouse head angle displacement w.r.t. intial  
                figure(h_angle);
                plot(d.time(selector, :), d.theta_filt_wrt_LED(selector, :), d.color_code{selector});
                hold on 
                vline(t_before_stim, 'b'); 
                hold off
                ylim([-10 180]); 
                yticks([0 90 180]); 
                xlim([0 -30 60]); 
                title([num2str(intensity), '; ', num2str(mouse_id)]);
                
                % Wait for user press any key
                pause(); 
               
            end 
        end
    end  
end 


%% Analysis of heading relative to LED, grouped by irradiance and experimental distribution
for i = 1:size(d, 1)
    [~, initial_heading_frame] = min(abs(d.time(i,:) - t_before_stim));
    [~, final_heading_frame] = min(abs(d.time(i,:) - t_after_stim));
    d.initial_heading(i) = d.theta_filt_wrt_LED(i, initial_heading_frame);
    d.final_heading(i) = d.theta_filt_wrt_LED(i, final_heading_frame);
end 

h_factor = figure;
h_angle = tiledlayout(length(unique_distributions),n_intensities, 'TileSpacing', 'compact');
h_angle.XLabel.String = 'Time (sec. from light onset)';
h_angle.YLabel.String = ['Head angle (', char(176), ')'];
h_angle.XLabel.FontSize = 10;
h_angle.YLabel.FontSize = 10;


h_theta = figure;
h_angle = tiledlayout(length(unique_distributions),n_intensities, 'TileSpacing', 'compact');
h_angle.XLabel.String = 'Time (sec. from light onset)';
h_angle.YLabel.String = ['Head angle (', char(176), ')'];
h_angle.XLabel.FontSize = 10;
h_angle.YLabel.FontSize = 10;

h_responders = figure;
d.turned = sum(d.theta_filt_wrt_LED >= 90, 2) > 0;
xspacing = 1:length(intensities);
jitter = [-0.1, 0.1]; 

for i = 1:length(unique_distributions)
    
    d_dis_orig = d(strcmp(d.distribution, unique_distributions{i}), :);
    d_dis = d_dis_orig;
    
    % Randomly select k mice from the distribution with unbroken tracking
    % data during test phase of experiment
    random_subset = false;
    if random_subset
        
        k = 1;
        
        for m = 1:size(d_dis, 1)
            trial_time = d_dis.time(m,:) < t_before_stim+60 & d_dis.time(m,:) > t_before_stim;
            d_dis.unbroken_track(m) = sum(~isnan(d_dis.time(m, trial_time)) ...
            ./ d_dis.framerate(m)) > 55;
        end
        cur_mouse_ids = unique(d_dis.mouse_id(d_dis.unbroken_track));

        rand_mouse_ids =  datasample(cur_mouse_ids, min([k, length(cur_mouse_ids)]),'Replace',false);
        
        mouse_id_rows = zeros(1, size(d_dis, 1));
        for m = 1:size(d_dis, 1)
            mouse_id_rows(m) = ~isempty(intersect(rand_mouse_ids, d_dis.mouse_id(m)));
        end

        d_dis = d_dis(logical(mouse_id_rows), :);
    end 
    
    % Set up mouse color code
    unique_mouse_ids = unique(d_dis.mouse_id);
    colorblind_colormap = load('colorblind_colormap.mat');
    colorblind_colormap = colorblind_colormap.colorblind;
    colorblind_colormap = [colorblind_colormap; [0 0 0]; [0 0 0]; [0 0 0]]; 
    colorblind_colormap = vertcat(colorblind_colormap(1:2, :), colorblind_colormap(4:end, :));
    mouse_color_code = colorblind_colormap(1:length(unique_mouse_ids), :);
    d_dis.mouse_color_code = mouse_color_code(grp2idx(d_dis.mouse_id), :);
    
    
    for j = 1:length(intensities)
        
        cur_d = d_dis(d_dis.intensity == intensities(j), :);
        cur_d_orig = d_dis_orig(d_dis_orig.intensity == intensities(j), :);
        
        % Factor plots (comparing initial and final heading)
        figure(h_factor);
        nexttile;
        makeFactorPlot([cur_d.initial_heading, cur_d.final_heading], ...
            h_factor, false, false, true, true, true, false, cur_d.color_code{1}, 'k', 3, [0.5, 0.5, 0.5]);
        if i == 1; title(num2str(intensities(j))); end
        ylim([-5 180]);
        yticks([0:90:180]);
        set(gca, 'FontSize', 8);
        xticklabels({'0', '60'});
    
        
        % Angle plots 
        figure(h_theta);
        nexttile;
        ylim([-5 180]);
        
        if random_subset
            
            for t = 1:size(cur_d_orig, 1)
                plot(cur_d_orig.time(t,:)-t_before_stim, cur_d_orig.theta_filt_wrt_LED(t,:),...
                'LineWidth', 0.25, 'Color', [0.5 0.5 0.5]);
                hold on 
            end
        end
        
        for t = 1:size(cur_d, 1)
            plot(cur_d.time(t,:)-t_before_stim, cur_d.theta_filt_wrt_LED(t,:),...
                'LineWidth', 0.5, 'Color', cur_d.mouse_color_code(t,:));
            hold on 
             
        end 
        
        ylim([-5 180]);
        xlim([0 90]);
        yticks(0:90:180); 
        xticks(0:30:90);
        set(gca, 'FontSize', 8, 'TickDir', 'out', 'box', 'off');
        
        if i == 1
            if intensities(j) == 0
                title('dark');
            else 
                title([num2str(intensities(j)), ' log units']);
            end
        end
        
        if j == 1; legend(num2str(unique(cur_d.mouse_id)), 'Location', 'northwest', 'FontSize', 6); end
        
        
        %Proportion responders (turned >= 90 degrees from LED) 
        figure(h_responders); 
        plot(xspacing(j)+jitter(i), mean(cur_d.turned*100), '.', 'Color', cur_d.color_code{1}, 'MarkerSize', 16);
        hold on 
        cis = bootci(1000, @mean, cur_d.turned*100);
        plot([xspacing(j)+jitter(i), xspacing(j)+jitter(i)], cis, 'Color', cur_d.color_code{1}, 'LineWidth', 1);

        set(h_responders, 'Position', [0   649   277   174]); 
        ylabel('Responders (%)');  
        xticks(xspacing);
        if length(xspacing) > 1
            xspacingdiff = xspacing(2) - xspacing(1); 
            xlim([min(xspacing) - 0.3 * xspacingdiff, max(xspacing) + 0.3 * xspacingdiff]);
        else
            xspacingdiff = 0;
        end 
        xticklabels(intensities);
        xlabel('Irradiance (log units)'); 
        set(gca, 'FontSize', 8, 'TickDir', 'out');
        box off

        
        
    end
end
set(h_factor, 'Position', [-929   309   368   250]);
set(h_theta, 'Renderer', 'painters', 'Position', [-1709         202         938         356]);
















