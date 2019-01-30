function analyzeData()
INVALID_EYE_COORD = -2^15;
MILLIS_PER_DRAW = 20;
STARTING_TRIAL = 1; 
LAST_FIXATION_BEFORE_SACCADE_DUR = 2000;
plot_title = [];
mutex = figure('visible', 'off');
DO_ANIMATE = true;

% These define a segmentation
CONDS = {'1'};
TRIAL_START_TRIGGER = '1';
PURSUIT_START_TRIGGER = '3';
PURSUIT_END_TRIGGERS = {'110', '111', '112', '120', '121', '122', '130', '131', '132', '140', '141', '142', '151', '152', '150'};
TRIAL_END_TRIGGER = '4';
    function trigger = getTrialTrigger(curr_trial_behavioral_data)
        trigger = '1';      
    end

% define data folders nemes
DATA_FOLDER = fullfile('Data');
EYE_DATA_PATH = fullfile(DATA_FOLDER, 'EyeData');
BEHAVIORAL_DATA_PATH = fullfile(DATA_FOLDER, 'BehavioralData');
ANALYSIS_STRUCT_PATH = fullfile(DATA_FOLDER);
eye_data_files_struct = dir(fullfile(EYE_DATA_PATH, 'Seye_example.mat')); % Eye data file converted to mat
expdata_files_struct = dir(fullfile(BEHAVIORAL_DATA_PATH, 'Sexpdata_example.mat')); % Behavioural mat file 
loaded_analysis_struct = load(fullfile(ANALYSIS_STRUCT_PATH, 's_analysis_struct_example.mat'));% Structure of detected fixations (based on Engbert and Kligel 2003) 
analysis_structs_cell_arr = loaded_analysis_struct.analysis_struct;
gui= figure('Visible', 'on', 'name', 'pursuit simulator', 'NumberTitle', 'off', 'units', 'pixels', ...
    'MenuBar', 'none', ...  
    'color', [0.8, 0.8, 0.8]);
main_win_axes = axes();
menu_bar= uimenu(gui,'Label','Action');
uimenu(menu_bar, 'Label', 'replay', 'callback', @setReplay);
uimenu(menu_bar, 'Label', 'skip', 'callback', @setSkip);
uimenu(menu_bar, 'Label', 'mark next fixation', 'callback', @markNextFixation);
uimenu(menu_bar, 'Label', 'mark previous fixation', 'callback', @markPrevFixation);
uimenu(menu_bar, 'Label', 'next trial', 'callback', @nextTrial);
analysis_struct_augmented = cell(1, numel(eye_data_files_struct));
for subject_i= 1:numel(eye_data_files_struct)
    disp(['Analyzing subject #', num2str(subject_i), ':']);
    disp('Loading eye data');
    % loading eye data file for current subject
    loaded_eye_tracking_data_mat_struct = load(fullfile(EYE_DATA_PATH, eye_data_files_struct(subject_i).name));
    %eye_tracking_data_mat = loaded_eye_tracking_data_mat_struct.eye_tracking_data_mat;
    loaded_eye_tracking_data_mat_struct.eye_tracking_data_structs
    eye_tracking_data_mat = loaded_eye_tracking_data_mat_struct.eye_tracking_data_structs{1, 1};
    disp('Loading expdata');
    % loading expdata file for current subject
    loaded_expdata_files_struct = load(fullfile(BEHAVIORAL_DATA_PATH, expdata_files_struct(subject_i).name));
    expdata = loaded_expdata_files_struct.EXPDATA;
        
    % eyelink messages struct iterator    
    pursuit_start_times = NaN(1, numel(expdata.trials));
    pursuit_end_times = NaN(1, numel(expdata.trials));        
    
    % ----------------------------
    % -----segmentizing data------
    % ----------------------------
    msg_i = 1;
    segmentized_eye_data = {};
    search_phase = 1;
    trial_idx = 1;
    disp('Segmentizing eye data');    
    while ~isempty(eye_tracking_data_mat.messages(msg_i).message)
        if search_phase == 1 && strcmp(eye_tracking_data_mat.messages(msg_i).message, TRIAL_START_TRIGGER)
            trial_start_time_i = find(eye_tracking_data_mat.gazeRight.time == eye_tracking_data_mat.messages(msg_i).time);
            search_phase = 2;
            continue;
        elseif search_phase == 2 && strcmp(eye_tracking_data_mat.messages(msg_i + 1).message, PURSUIT_START_TRIGGER)
            pursuit_start_times(trial_idx) = find(eye_tracking_data_mat.gazeRight.time == eye_tracking_data_mat.messages(msg_i + 1).time) - trial_start_time_i + 1;
            search_phase = 3;
        elseif search_phase == 3 && any(cellfun(@(s) strcmp(eye_tracking_data_mat.messages(msg_i + 1).message, s), PURSUIT_END_TRIGGERS))
            pursuit_end_times(trial_idx) = find(eye_tracking_data_mat.gazeRight.time == eye_tracking_data_mat.messages(msg_i + 1).time) - trial_start_time_i + 1;
            search_phase = 4;
        elseif search_phase == 4 && (any(cellfun(@(s) strcmp(eye_tracking_data_mat.messages(msg_i + 1).message, s), CONDS)) || ...
                strcmp(eye_tracking_data_mat.messages(msg_i + 1).message, TRIAL_END_TRIGGER))
            trial_end_time_i = find(eye_tracking_data_mat.gazeRight.time == (eye_tracking_data_mat.messages(msg_i + 1).time));
            if isempty(trial_start_time_i) || isempty(trial_end_time_i)
                segmentized_eye_data{trial_idx}.gazeRight = [];
                segmentized_eye_data{trial_idx}.gazeLeft = [];
            else
                trial_gaze_data_is = trial_start_time_i:trial_end_time_i;
                segmentized_eye_data{trial_idx}.gazeRight = [eye_tracking_data_mat.gazeRight.x(trial_gaze_data_is);
                    eye_tracking_data_mat.gazeRight.y(trial_gaze_data_is)]';
                segmentized_eye_data{trial_idx}.gazeRight(segmentized_eye_data{trial_idx}.gazeRight == INVALID_EYE_COORD) = NaN;
                segmentized_eye_data{trial_idx}.gazeLeft = [eye_tracking_data_mat.gazeLeft.x(trial_gaze_data_is);
                    eye_tracking_data_mat.gazeLeft.y(trial_gaze_data_is)]';
                segmentized_eye_data{trial_idx}.gazeLeft(segmentized_eye_data{trial_idx}.gazeRight == INVALID_EYE_COORD) = NaN;
            end
            search_phase = 1;           
            trial_idx = trial_idx + 1;
        elseif strcmp(eye_tracking_data_mat.messages(msg_i + 1).message, '200') || ...               
                strcmp(eye_tracking_data_mat.messages(msg_i + 1).message, '68')
            search_phase = 1;
        end
        
        msg_i = msg_i + 1;        
    end            
        
    analysis_struct_iterators = ones(1, numel(CONDS));
    trials_it = 1;
    for triSTARTING_TRIALal_i = 1: - 1       
        if isempty(expdata.trials(trial_i).trial_duration) || strcmp(expdata.trials(trial_i).subject_response, 'pursuit brake')
            continue;
        end
        cond_i = find(cellfun(@(cond) strcmp(getTrialTrigger(expdata.trials(trial_i)), cond), CONDS));        
        analysis_struct_iterators(cond_i) = analysis_struct_iterators(cond_i) + 1;  
        trials_it = trials_it + 1;        
    end
    analysis_struct_augmented{subject_i} = analysis_structs_cell_arr{subject_i}.c1; % TEMP
        
    for trial_i= STARTING_TRIAL:numel(expdata.trials)

        curr_trial_behavioral_data = expdata.trials(trial_i);
        % relevant for subject 102/202
        if isempty(curr_trial_behavioral_data.trial_duration) || strcmp(curr_trial_behavioral_data.subject_response, 'pursuit brake')
            continue;
        end
        
        curr_trial_trigger = getTrialTrigger(curr_trial_behavioral_data);
        cond_i = find(cellfun(@(cond) strcmp(curr_trial_trigger, cond), CONDS));

        % -----------------------------
        % -----drawing next trial------
        % -----------------------------
        curr_trial_fixations_analysis_struct = analysis_structs_cell_arr{subject_i}.(['c',curr_trial_trigger]).fixations(analysis_struct_iterators(cond_i));        
        
        % extracting trial's parameters <<FROM EXPDATA>>
        drawn_pursuit_dur = expdata.info.experiment_parameters.travel_durations(abs(expdata.info.experiment_parameters.travel_durations - curr_trial_behavioral_data.pursuit_duration) < 0.005);
        
        ppd = 60; % pixels per degree
        if strcmp(curr_trial_behavioral_data.initial_target_side, 'left')
            target_initial_x = expdata.info.experiment_parameters.target_initial_distance_from_edge * ppd;
            target_final_x = target_initial_x + expdata.info.experiment_parameters.target_speed*ppd*drawn_pursuit_dur;
            if strcmp(curr_trial_behavioral_data.retinotopic_placeholder_side, 'up')
                upper_placeholder_x = target_final_x + curr_trial_behavioral_data.target_to_retino_vec_radius * ppd * cos(deg2rad(curr_trial_behavioral_data.target_to_retino_vec_angle));
                upper_placeholder_y = expdata.info.lab_setup.screen_height_in_pixels / 2 - curr_trial_behavioral_data.target_to_retino_vec_radius * ppd * sin(deg2rad(curr_trial_behavioral_data.target_to_retino_vec_angle));
                lower_placeholder_x = target_final_x + curr_trial_behavioral_data.target_to_spatio_vec_radius * ppd * cos(deg2rad(curr_trial_behavioral_data.target_to_spatio_vec_angle));
                lower_placeholder_y = expdata.info.lab_setup.screen_height_in_pixels / 2 + curr_trial_behavioral_data.target_to_spatio_vec_radius * ppd * sin(deg2rad(curr_trial_behavioral_data.target_to_spatio_vec_angle));
                upper_placeholder_color = expdata.info.experiment_parameters.retinotopic_placeholder_color;
                lower_placeholder_color = expdata.info.experiment_parameters.spaciotopic_placeholder_color;
            else 
                upper_placeholder_x = target_final_x + curr_trial_behavioral_data.target_to_spatio_vec_radius * ppd * cos(deg2rad(curr_trial_behavioral_data.target_to_spatio_vec_angle));
                upper_placeholder_y = expdata.info.lab_setup.screen_height_in_pixels / 2 - curr_trial_behavioral_data.target_to_spatio_vec_radius * ppd * sin(deg2rad(curr_trial_behavioral_data.target_to_spatio_vec_angle));
                lower_placeholder_x = target_final_x + curr_trial_behavioral_data.target_to_retino_vec_radius * ppd * cos(deg2rad(curr_trial_behavioral_data.target_to_retino_vec_angle));
                lower_placeholder_y = expdata.info.lab_setup.screen_height_in_pixels / 2 + curr_trial_behavioral_data.target_to_retino_vec_radius * ppd * sin(deg2rad(curr_trial_behavioral_data.target_to_retino_vec_angle));
                upper_placeholder_color = expdata.info.experiment_parameters.spaciotopic_placeholder_color;
                lower_placeholder_color = expdata.info.experiment_parameters.retinotopic_placeholder_color;
            end

        elseif strcmp(curr_trial_behavioral_data.initial_target_side, 'right')
            target_initial_x = expdata.info.lab_setup.screen_width_in_pixels - expdata.info.experiment_parameters.target_initial_distance_from_edge * ppd;
            target_final_x = target_initial_x - expdata.info.experiment_parameters.target_speed*ppd*drawn_pursuit_dur;
            if strcmp(curr_trial_behavioral_data.retinotopic_placeholder_side, 'up')
                upper_placeholder_x = target_final_x - curr_trial_behavioral_data.target_to_retino_vec_radius * ppd * cos(deg2rad(curr_trial_behavioral_data.target_to_retino_vec_angle));
                upper_placeholder_y = expdata.info.lab_setup.screen_height_in_pixels / 2 - curr_trial_behavioral_data.target_to_retino_vec_radius * ppd * sin(deg2rad(curr_trial_behavioral_data.target_to_retino_vec_angle));
                lower_placeholder_x = target_final_x - curr_trial_behavioral_data.target_to_spatio_vec_radius * ppd * cos(deg2rad(curr_trial_behavioral_data.target_to_spatio_vec_angle));
                lower_placeholder_y = expdata.info.lab_setup.screen_height_in_pixels / 2 + curr_trial_behavioral_data.target_to_spatio_vec_radius * ppd * sin(deg2rad(curr_trial_behavioral_data.target_to_spatio_vec_angle));
                upper_placeholder_color = expdata.info.experiment_parameters.retinotopic_placeholder_color;
                lower_placeholder_color = expdata.info.experiment_parameters.spaciotopic_placeholder_color;
            else
                upper_placeholder_x = target_final_x - curr_trial_behavioral_data.target_to_spatio_vec_radius * ppd * cos(deg2rad(curr_trial_behavioral_data.target_to_spatio_vec_angle));
                upper_placeholder_y = expdata.info.lab_setup.screen_height_in_pixels / 2 - curr_trial_behavioral_data.target_to_spatio_vec_radius * ppd * sin(deg2rad(curr_trial_behavioral_data.target_to_spatio_vec_angle));
                lower_placeholder_x = target_final_x - curr_trial_behavioral_data.target_to_retino_vec_radius * ppd * cos(deg2rad(curr_trial_behavioral_data.target_to_retino_vec_angle));
                lower_placeholder_y = expdata.info.lab_setup.screen_height_in_pixels / 2 + curr_trial_behavioral_data.target_to_retino_vec_radius * ppd * sin(deg2rad(curr_trial_behavioral_data.target_to_retino_vec_angle));
                upper_placeholder_color = expdata.info.experiment_parameters.spaciotopic_placeholder_color;
                lower_placeholder_color = expdata.info.experiment_parameters.retinotopic_placeholder_color;
            end
        end
      
        IS_TRIAL_FINISHED = false;
        DO_REPLAY = false;
        DO_SKIP = false;
        while 1
            if DO_ANIMATE
                axes(main_win_axes);
                cla;
                set(gca,'yDir', 'reverse', 'xAxisLocation', 'top');
                if ~isempty(plot_title)
                    delete(plot_title);
                end
                if curr_trial_behavioral_data.is_gabor_tilted
                    plot_title = text(0,0,['trial #', num2str(trial_i), '; gabor appearance side: ', curr_trial_behavioral_data.gabor_appearance_side, '; block type: ', expdata.blocks(curr_trial_behavioral_data.block_number).block_type, '; retinotopic placeholder side: ', curr_trial_behavioral_data.retinotopic_placeholder_side, '; has subject responded: ', num2str(strcmp(curr_trial_behavioral_data.subject_response, 'saccade'))]);
                else
                    plot_title =  text(0,0,['trial #', num2str(trial_i), '; gabor appearance side: nowhere; block type: ', expdata.blocks(curr_trial_behavioral_data.block_number).block_type, '; retinotopic placeholder side: ', curr_trial_behavioral_data.retinotopic_placeholder_side, '; has subject responded: ', num2str(strcmp(curr_trial_behavioral_data.subject_response, 'saccade'))]);
                end
                set(gca, 'title', plot_title);

                hold on;

                % draw target at it's initial position
                plot(target_initial_x, expdata.info.lab_setup.screen_height_in_pixels / 2, 'ok', 'markersize', 20);
                % draw the upper placeholder at it's final position
                plot(upper_placeholder_x, upper_placeholder_y, 'o', 'markersize', 20, 'color', upper_placeholder_color);
                % draw the lower placeholder at it's final position
                plot(lower_placeholder_x, lower_placeholder_y, 'o', 'markersize', 20, 'color', lower_placeholder_color);

                % calculate the target's travel distance per frame
                if strcmp(curr_trial_behavioral_data.initial_target_side, 'left')
                    target_travel_dist_in_pixels_per_draw = expdata.info.experiment_parameters.target_speed * ppd / 1000 * MILLIS_PER_DRAW;
                else
                    target_travel_dist_in_pixels_per_draw = -expdata.info.experiment_parameters.target_speed * ppd / 1000 * MILLIS_PER_DRAW;
                end
                time_until_pursuit_onset = (curr_trial_behavioral_data.wait_for_fixation_duration + ...
                    curr_trial_behavioral_data.initial_fixation_on_target_duration + ...
                    curr_trial_behavioral_data.cue_duration + ...
                    curr_trial_behavioral_data.interval_duration_between_cue_and_pursuit)*1000;
                prev_target_x = target_initial_x;

                disp('Simulating trial');
                for sample_i = 1:MILLIS_PER_DRAW:length(segmentized_eye_data{trials_it}.gazeLeft)
                    if sample_i > time_until_pursuit_onset
                        plot(prev_target_x, expdata.info.lab_setup.screen_height_in_pixels / 2, 'ow', 'markersize', 20);
                        plot(prev_target_x + target_travel_dist_in_pixels_per_draw, expdata.info.lab_setup.screen_height_in_pixels / 2, 'ok', 'markersize', 20);
                        prev_target_x = prev_target_x + target_travel_dist_in_pixels_per_draw;
                    end

                    % draw subject's eye positions
                    if ~isnan(segmentized_eye_data{trials_it}.gazeLeft(sample_i,1))
                        plot(segmentized_eye_data{trials_it}.gazeLeft(sample_i,1), segmentized_eye_data{trials_it}.gazeLeft(sample_i,2), '.g', 'markersize', 5);
                    end
                    if ~isnan(segmentized_eye_data{trials_it}.gazeRight(sample_i,1))
                        plot(segmentized_eye_data{trials_it}.gazeRight(sample_i,1), segmentized_eye_data{trials_it}.gazeRight(sample_i,2), '.m', 'markersize', 5);
                    end
                    pause(MILLIS_PER_DRAW/1000);
                    if DO_SKIP || DO_REPLAY
                        break;
                    end
                end

                if DO_SKIP
                    DO_SKIP = false;
                    analysis_struct_iterators(cond_i) = analysis_struct_iterators(cond_i) + 1;
                    trials_it = trials_it + 1;
                    break;
                elseif DO_REPLAY
                    DO_REPLAY = false;
                    continue;
                end            
            end
        
            fixations_nr = numel(curr_trial_fixations_analysis_struct.fixations_durations);
            last_fixation_before_pursuit_idx = fixations_nr - find(curr_trial_fixations_analysis_struct.fixations_durations(end:-1:1) > LAST_FIXATION_BEFORE_SACCADE_DUR, 1) + 1;
            if isempty(last_fixation_before_pursuit_idx) || last_fixation_before_pursuit_idx == fixations_nr
                response_fixation_idx = NaN;
            else
                response_fixation_idx = last_fixation_before_pursuit_idx + 1;
            end
            
            if DO_ANIMATE
                fixations_plots = NaN(fixations_nr, 2);            
                for fixation_i = 1:fixations_nr
                    if abs(curr_trial_fixations_analysis_struct.fixations_coordinates_left(fixation_i,1) - INVALID_EYE_COORD) > 1e-5 && ...
                            abs(curr_trial_fixations_analysis_struct.fixations_coordinates_right(fixation_i,1) - INVALID_EYE_COORD) > 1e-5
                        mean_fixation_coords = mean([curr_trial_fixations_analysis_struct.fixations_coordinates_left(fixation_i,:);
                            curr_trial_fixations_analysis_struct.fixations_coordinates_right(fixation_i,:)]);
                        if fixation_i == response_fixation_idx
                            left_eye_fixation_plot_desc = 'ob';
                            right_eye_fixation_plot_desc = 'ob';
                        else
                            left_eye_fixation_plot_desc = 'og';
                            right_eye_fixation_plot_desc = 'om';
                        end
                        fixations_plots(fixation_i, 1) = plot(gca, curr_trial_fixations_analysis_struct.fixations_coordinates_left(fixation_i,1), curr_trial_fixations_analysis_struct.fixations_coordinates_left(fixation_i,2), left_eye_fixation_plot_desc , 'markersize', 10);
                        fixations_plots(fixation_i, 2) = plot(gca, curr_trial_fixations_analysis_struct.fixations_coordinates_right(fixation_i,1), curr_trial_fixations_analysis_struct.fixations_coordinates_right(fixation_i,2), right_eye_fixation_plot_desc , 'markersize', 10);
                        enum_valid_eyes = 0;
                    elseif abs(curr_trial_fixations_analysis_struct.fixations_coordinates_left(fixation_i,1) - INVALID_EYE_COORD) > 1e-5
                        mean_fixation_coords = curr_trial_fixations_analysis_struct.fixations_coordinates_left(fixation_i,:);
                        if fixation_i == response_fixation_idx
                            fixation_plot_desc = 'ob';
                        else
                            fixation_plot_desc = 'og';
                        end
                        fixations_plots(fixation_i, 1) = plot(gca, curr_trial_fixations_analysis_struct.fixations_coordinates_left(fixation_i,1), curr_trial_fixations_analysis_struct.fixations_coordinates_left(fixation_i,2), fixation_plot_desc , 'markersize', 10);
                        enum_valid_eyes = 1;
                    elseif abs(curr_trial_fixations_analysis_struct.fixations_coordinates_right(fixation_i,1) - INVALID_EYE_COORD) > 1e-5
                        mean_fixation_coords = curr_trial_fixations_analysis_struct.fixations_coordinates_right(fixation_i,:);
                        if fixation_i == response_fixation_idx
                            fixation_plot_desc = 'ob';
                        else
                            fixation_plot_desc = 'ob';
                        end
                        fixations_plots(fixation_i, 2) = plot(gca, curr_trial_fixations_analysis_struct.fixations_coordinates_right(fixation_i,1), curr_trial_fixations_analysis_struct.fixations_coordinates_right(fixation_i,2), fixation_plot_desc, 'markersize', 10);
                        enum_valid_eyes = 2;
                    else
                        continue;
                    end
                    text(mean_fixation_coords(1)-15, mean_fixation_coords(2), num2str(fixation_i), 'color', [1 0 0]);
                    text(mean_fixation_coords(1)+25, mean_fixation_coords(2), num2str(curr_trial_fixations_analysis_struct.fixations_durations(fixation_i)), 'color', [0 0 0]);
                end

                IS_TRIAL_FINISHED = true;       
                disp('You are eyeballing');
                waitfor(mutex);      
            end
            
            if DO_REPLAY
                DO_REPLAY = false;               
            else                                                
                analysis_struct_augmented{subject_i}.fixations(analysis_struct_iterators(cond_i)).block_type = expdata.blocks(curr_trial_behavioral_data.block_number).block_type;
                analysis_struct_augmented{subject_i}.fixations(analysis_struct_iterators(cond_i)).response_fixation_index = response_fixation_idx;
                analysis_struct_augmented{subject_i}.fixations(analysis_struct_iterators(cond_i)).was_gabor_tilted = curr_trial_behavioral_data.is_gabor_tilted;
                analysis_struct_augmented{subject_i}.fixations(analysis_struct_iterators(cond_i)).upper_place_holder_location = [upper_placeholder_x, upper_placeholder_y];
                analysis_struct_augmented{subject_i}.fixations(analysis_struct_iterators(cond_i)).lower_place_holder_location = [lower_placeholder_x, lower_placeholder_y];
                if curr_trial_behavioral_data.is_gabor_tilted
                    analysis_struct_augmented{subject_i}.fixations(analysis_struct_iterators(cond_i)).tilted_gabor_side = curr_trial_behavioral_data.gabor_appearance_side;
                else
                    analysis_struct_augmented{subject_i}.fixations(analysis_struct_iterators(cond_i)).tilted_gabor_side = 'nowhere';
                end
                analysis_struct_augmented{subject_i}.saccades(analysis_struct_iterators(cond_i)).catchup_saccades_indices = ...
                    find(analysis_structs_cell_arr{subject_i}.(['c',curr_trial_trigger]).saccades(analysis_struct_iterators(cond_i)).onsets < pursuit_end_times(trials_it) & pursuit_start_times(trials_it) < analysis_structs_cell_arr{subject_i}.(['c',curr_trial_trigger]).saccades(analysis_struct_iterators(cond_i)).offsets);
                analysis_struct_augmented{subject_i}.raw_data(analysis_struct_iterators(cond_i)).pursuit_start_time = pursuit_start_times(trials_it);
                analysis_struct_augmented{subject_i}.raw_data(analysis_struct_iterators(cond_i)).pursuit_end_time = pursuit_end_times(trials_it);
                analysis_struct_augmented{subject_i}.raw_data(analysis_struct_iterators(cond_i)).left_eye.x = segmentized_eye_data{trials_it}.gazeLeft(:,1);
                analysis_struct_augmented{subject_i}.raw_data(analysis_struct_iterators(cond_i)).left_eye.y = segmentized_eye_data{trials_it}.gazeLeft(:,2);
                analysis_struct_augmented{subject_i}.raw_data(analysis_struct_iterators(cond_i)).right_eye.x = segmentized_eye_data{trials_it}.gazeRight(:,1);
                analysis_struct_augmented{subject_i}.raw_data(analysis_struct_iterators(cond_i)).right_eye.y = segmentized_eye_data{trials_it}.gazeRight(:,2);
                analysis_struct_augmented{subject_i}.raw_data(analysis_struct_iterators(cond_i)).pursuit_dot_start_location = [target_initial_x, expdata.info.lab_setup.screen_height_in_pixels / 2];
                analysis_struct_augmented{subject_i}.raw_data(analysis_struct_iterators(cond_i)).pursuit_dot_final_location = [target_final_x, expdata.info.lab_setup.screen_height_in_pixels / 2];
                analysis_struct_iterators(cond_i) = analysis_struct_iterators(cond_i) + 1;
                trials_it = trials_it + 1;
                if DO_ANIMATE
                    disp('Saving intermediate analysis data');
                    %save(fullfile(DATA_FOLDER, 'analysis_struct_augmented.mat'), 'analysis_struct_augmented');
                end
                break;
            end
        end
    end
    
    close all
    disp(['done with subject #', num2str(subject_i), '.']);    
end

% to save the data
%disp('Saving final analysis data');
%save(fullfile(DATA_FOLDER, 'analysis_struct_augmented.mat'), 'analysis_struct_augmented');
%disp('Done.');


    function setReplay(~,~)
        DO_REPLAY = true; 
        if IS_TRIAL_FINISHED
            delete(mutex);  
            mutex = figure('visible', 'off');            
        end
    end

    function setSkip(~,~)
        DO_SKIP = true;
    end

    function markNextFixation(~,~)        
        if ~IS_TRIAL_FINISHED || isnan(response_fixation_idx)
            return;                    
        else        
            if enum_valid_eyes == 1
                if response_fixation_idx > 0
                    set(fixations_plots(response_fixation_idx, 1), 'Color', [0, 1, 0]);
                end
                response_fixation_idx = response_fixation_idx + 1;
                if response_fixation_idx > fixations_nr
                    response_fixation_idx = NaN;
                else
                    set(fixations_plots(response_fixation_idx, 1), 'color', [0, 0, 1]);
                end
            elseif enum_valid_eyes == 2    
                if response_fixation_idx > 0
                    set(fixations_plots(response_fixation_idx, 2), 'color', [1, 0, 1]);
                end
                response_fixation_idx = response_fixation_idx + 1;
                if response_fixation_idx > fixations_nr
                    response_fixation_idx = NaN;
                else
                    set(fixations_plots(response_fixation_idx, 1), 'color', [0, 0, 1]);
                end
            else
                if response_fixation_idx > 0
                    set(fixations_plots(response_fixation_idx, 1), 'color', [0, 1, 0]);
                    set(fixations_plots(response_fixation_idx, 2), 'color', [1, 0, 1]);
                end
                response_fixation_idx = response_fixation_idx + 1;
                if response_fixation_idx > fixations_nr
                    response_fixation_idx = NaN;
                else
                    set(fixations_plots(response_fixation_idx, 1), 'color', [0, 0, 1]);
                    set(fixations_plots(response_fixation_idx, 2), 'color', [0, 0, 1]);
                end
            end
        end
    end

    function markPrevFixation(~,~)        
        if ~IS_TRIAL_FINISHED || response_fixation_idx == 0
            return;                    
        else        
            if enum_valid_eyes == 1
                if ~isnan(response_fixation_idx)
                    set(fixations_plots(response_fixation_idx, 1), 'color', [0, 1, 0]);
                else
                end
                response_fixation_idx = response_fixation_idx - 1;
                if response_fixation_idx > 0                    
                    set(fixations_plots(response_fixation_idx, 1), 'color', [0, 0, 1]);
                end
            elseif enum_valid_eyes == 2     
                if ~isnan(response_fixation_idx)
                    set(fixations_plots(response_fixation_idx, 2), 'color', [1, 0, 1]);
                end
                response_fixation_idx = response_fixation_idx - 1;
                if response_fixation_idx > 0                    
                    set(fixations_plots(response_fixation_idx, 1), 'color', [0, 0, 1]);
                end
            else
                if ~isnan(response_fixation_idx)
                    set(fixations_plots(response_fixation_idx, 1), 'color', [0, 1, 0]);
                    set(fixations_plots(response_fixation_idx, 2), 'color', [1, 0, 1]);
                    response_fixation_idx = response_fixation_idx - 1;
                else
                    response_fixation_idx = fixations_nr;
                end                
                if response_fixation_idx > 0                  
                    set(fixations_plots(response_fixation_idx, 1), 'color', [0, 0, 1]);
                    set(fixations_plots(response_fixation_idx, 2), 'color', [0, 0, 1]);
                end
            end
        end
    end

    function nextTrial(~,~)
        if IS_TRIAL_FINISHED
            delete(mutex);  
            mutex = figure('visible', 'off');            
        end
    end
end


