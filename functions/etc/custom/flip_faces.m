%%  SCRIPT TO AUTOMATICALLY/MANUALLY FLIP FACES TIMESERIES

% Outputs:
%    flip.mat       a file flip.mat, specifying whether subjects should be 
%                   flipped (1 = flip, 0 = don't flip)
%
%    occ_ts.png     plots of the occipital regions used to calculate the 
%                   flip value,stored in the analysis subfolder

close all
clear
clc

% Run the study configuration
p = opm_study_config_faces();

% Overwrite?
clobber = 0;

% Define idx to look at
ffg = 56;

% Specify mode
mode = 'manual';

% Output directory
out_dir = [p.directories.base_dir '/analysis/' p.pipeline_name '/face_flipping'];
filename = [out_dir '/' 'flip_l.mat'];

% Define subjects to look at
subjects = 1:size(p.subject_data, 1);
num_subjects = length(subjects);

% Get time vector
time = get_time(p);
num_timepoints = length(time);

% Preallocate memory
flip = zeros(num_subjects, 1);
ts_ffg = zeros(1, num_timepoints, num_subjects);

% Loop over subjects
for i = 1:num_subjects

    % Index for p
    ss = subjects(i);
    fprintf('Loading data for subject %s %s %s\n', p.subject(ss), p.session(ss), p.run(ss))

    % Get ffg response
    [ts_ffg(:, :, i), no_trials] = get_timeseries(p, ffg, ss);
    if anynan(ts_ffg(:, :, i))
        % This means something went wrong during source modeling
        disp('Timeseries contains NaN value - breaking')
        flip(i) = NaN;
    else
        % Flip if needed 
        fprintf('Flipping timeseries for subject %s %s %s\n', p.subject(ss), p.session(ss), p.run(ss))
        
        if strcmp(mode, 'auto')
            flip(i) = flip_face(p, ss, ts_ffg(:, :, i));
        elseif strcmp(mode, 'manual')
            figure
            plot(time, ts_ffg(:,:, i))
            xline([0.1 0.17 0.25])
            xlim([-0.1 0.5])
            
            prompt = {'Needs flipping? y/n: '};
            dlgtitle = 'Input';
            fieldsize = [1 45];
            definput = {'n'};
            opts = []; 
            opts.WindowStyle = 'normal';
            answer = inputdlg(prompt, dlgtitle, fieldsize, definput, opts);
            close all

            % Check the user's response
            if ~isempty(answer)  % Ensure the user didn't close the dialog
                user_input = lower(answer{1});  % Convert to lowercase for case-insensitivity
                if strcmp(user_input, 'y')
                    disp('Flipping...');
                    flip(i) = 1;
                elseif strcmp(user_input, 'n')
                    flip(i) = 0;
                    disp('No flipping needed.');
                end
            end
        end
    end
end


% Save flip vector
%save(filename, 'flip')
