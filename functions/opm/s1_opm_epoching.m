function s1_opm_epoching(p, sub_num, clobber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s1_opm_epoching(p, sub_num, clobber)
%
% Script to read in OPM MEG task data, convert to BIDS format that's
% acceptable for FieldTrip, and epoch the data by task conditions. Noise
% data is also read in and converted, if it's been specified.
%
% Inputs:
%   p                   (structure) Study configuration structure
%   sub_num             (ind)       Subject number to process
%   clobber             (logical)   Whether or not to overwrite existing
%                                   files
%
% Outputs:
%   "task_name"         (directory) BIDS directory of task and noise data
%   "trial_cfg.mat"     (file)      FieldTrip trial configuration
%
% Created by:   Marlee M Vandewouw
% Created on:   January 20th, 2023
% Modified on:  August 14th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION

% Make derivatives directory
if ~exist(p.directories.sub_dir(p.subject(sub_num), p.session(sub_num)), 'dir')
    cmd = ['mkdir -p ' p.directories.sub_dir(p.subject(sub_num), p.session(sub_num))];
    system(cmd);
end

%% CONVERT RAW DATA TO BIDS FORMAT

% Only convert if doesn't exist, otherwise, link
if p.directory.bids
    
    %% LINK TASK DATA
    
    % Only re-link if the link doesn't exist, or if we need to re-link
    if ~exist(p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1
        
        % Link
        fprintf('   Linking task data...');
        
        % Unlink if necessary
        if exist(p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file')
            
            cmd = ['unlink ' p.directories.raw_pre(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)) '*'];
            system(cmd);
            
        end
        
        % Make the base directory
        if ~exist(p.directories.raw_dir(p.subject(sub_num), p.session(sub_num)), 'dir')
            cmd = ['mkdir -p ' p.directories.raw_dir(p.subject(sub_num), p.session(sub_num))];
            system(cmd);
        end        
        
        % Link all the files corresponding to the task
        cmd = ['ln -s ' p.directory.bids_dir '/' p.subject(sub_num) '/' p.session(sub_num) '/opm/*' p.task(sub_num) '*' p.run(sub_num) '* ' p.directories.raw_dir(p.subject(sub_num), p.session(sub_num)) '/'];
        system(cmd);  
        
        % Link the noise
        cmd = ['ln -s ' p.directory.bids_dir '/' p.subject(sub_num) '/' p.session(sub_num) '/opm/*noise*' p.noise_run(sub_num) '* ' p.directories.raw_dir(p.subject(sub_num), p.session(sub_num)) '/'];
        system(cmd);
        
        fprintf('done\n');
        
    end
    
    % Link the anatomical data, if it exists
    if exist([p.directory.bids_dir '/' p.subject(sub_num) '/' p.session(sub_num) '/anat'], 'dir')
        
        if ~exist(p.directories.anat_dir(sub_num), 'dir') || clobber == 1
            
            % Delete if necessary
            if exist(p.directories.anat_dir(sub_num), 'dir') 
                cmd = ['rm -r ' p.directories.anat_dir(sub_num)];
                system(cmd);
            end
            
            % Make the directory
            cmd = ['mkdir ' p.directories.anat_dir(sub_num)];
            system(cmd);
            
            % Link
            cmd = ['ln -s ' p.directory.bids_dir '/' p.subject(sub_num) '/' p.session(sub_num) '/anat/* ' p.directories.anat_dir(sub_num) '/'];
            system(cmd);
            
        end
        
    end
        

else
    
    %% CONVERT TASK DATA
    
    % Only convert it we need to
    if ~exist(p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1
        
        fprintf('   Converting task data...');
        
        % If exists, delete
        if exist(p.directories.raw_dir(p.subject(sub_num), p.session(sub_num)), 'dir')
            cmd = ['rm -r ' p.directories.raw_dir(p.subject(sub_num), p.session(sub_num)) '/*'];
            system(cmd);
        end

        % Check what format the cMEG data by checking the date
        scan_date = strsplit(p.subject_data.opm_data{sub_num}, '/');
        scan_date = strsplit(scan_date{end}, '_');
        scan_date = datetime(str2double(scan_date{1}),'ConvertFrom','yyyymmdd');

        % Prepare conversion structure
        info = [];
        info.filename = p.subject_data.opm_data{sub_num};
        info.pos_dir = p.directories.helmet_dir;
        info.output_dir = p.directories.raw_dir(p.subject(sub_num), p.session(sub_num));

        % Check which conversion to use
        if scan_date < datetime(20230208,'ConvertFrom','yyyymmdd')

            % Use original function
            convert_cMEG_data(info);

        else

            % Use new function
            convert_cMEG_data_tsv(info);

        end

        % Move files to their appropriate names
        movefile([p.directories.raw_dir(p.subject(sub_num), p.session(sub_num)) 'meg.bin'], p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        movefile([p.directories.raw_dir(p.subject(sub_num), p.session(sub_num)) 'meg.json'], p.raw_data.json(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        movefile([p.directories.raw_dir(p.subject(sub_num), p.session(sub_num)) 'positions.tsv'], p.raw_data.positions(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        movefile([p.directories.raw_dir(p.subject(sub_num), p.session(sub_num)) 'channels.tsv'], p.raw_data.channels(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        
        fprintf('done\n');
        
    end
    
    %% CONVERT HEAD MOTION DATA

    % Only run if doesn't exist or if overwriting
    if ~exist(p.headmotion.file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1

        fprintf('   Converting headmotion data...');
        
        % Check to see if file or folder
        if exist(p.subject_data.opm_data{sub_num}, 'dir')
            % Directory, so search for csv
            d = dir([p.subject_data.opm_data{sub_num} '/*csv']);
            csv = [d.folder '/' d.name];
        else
            d = strsplit(p.subject_data.opm_data{sub_num}, '/');
            d = strjoin(d(1:end-1), '/');
            d = dir([d '/*csv']);
        end

        % Check if empty
        if isempty(d)
            warning('No optitrack data found!');
        else

            % Get .csv file for optitrack data
            if length(d) > 1
                error('More than one optitrack csv file found: please see Marlee');
            else
                csv = [d.folder '/' d.name];
            end

            % Load in the MEG data
            cfg_loadds = [];
            cfg_loadds.dataset = p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num));
            cfg_loadds.headerformat = 'opm_fil_mjt';
            cfg_loadds.dataformat = 'opm_fil_mjt';
            cfg_loadds.eventformat = 'opm_fil_mjt';
            data = ft_preprocessing(cfg_loadds);

            % Parse optitrack using the FIL-OPMEG toolbox
            try
                opti_data = csv2mat_sm(csv);

                % Get head motion data (note, this script was minorly modified to
                % conform with SK data)
                try
                    [MovementDataOut, OPMdataOut] = syncOptitrackAndOPMdata(opti_data, data, 'TriggerChannelName', 'Trigger 2 [Z]');

                    % Instead of cutting OPM data, padd the HeadMotion data
                    num_samples_pad = OPMdataOut.sampleinfo(1) - 1;
                    HeadMotionData = MovementDataOut;
                    HeadMotionData.rigidbodies.data = [zeros(num_samples_pad, size(MovementDataOut.rigidbodies.data, 2)); MovementDataOut.rigidbodies.data];
                    f = fieldnames(MovementDataOut.markers);
                    for ff = 1:length(f)
                        for ii = 1:length(MovementDataOut.markers.(sprintf('%s', f{ff})))
                            HeadMotionData.markers.(sprintf('%s', f{ff}))(ii).data = [zeros(num_samples_pad, size(HeadMotionData.markers.(sprintf('%s', f{ff}))(ii).data, 2)); HeadMotionData.markers.(sprintf('%s', f{ff}))(ii).data];
                        end
                    end
                    HeadMotionData.time = data.time{1};

                    % Save
                    parsave_struct(p.headmotion.file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), HeadMotionData, 'HeadMotionData');
                catch
                    warning('Head motion data unable to be processed for subject %s, session %s, task %s, run %s\n', p.subject(sub_num), p.session(sub_num), p.task(sub_num), p.run(sub_num));
                end

            catch

                warning('Head motion data unable to be processed for subject %s, session %s, task %s, run %s\n', p.subject(sub_num), p.session(sub_num), p.task(sub_num), p.run(sub_num));

            end

        end
        
        fprintf('done\n');

    end
    
    %% CONVERT NOISE DATA
    
    % Only convert it we need to
    if ~exist(p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), 'noise'), 'file') || clobber == 1
        
        fprintf('   Converting noise data...');

        % Check what format the cMEG data by checking the date
        scan_date = strsplit(p.subject_data.opm_noise{sub_num}, '/');
        scan_date = strsplit(scan_date{end}, '_');
        scan_date = datetime(str2double(scan_date{1}),'ConvertFrom','yyyymmdd');

        % Prepare conversion structure
        info = [];
        info.filename = p.subject_data.opm_noise{sub_num};
        info.pos_dir = p.directories.helmet_dir;
        info.output_dir = p.directories.raw_dir(p.subject(sub_num), p.session(sub_num));

        % Check which conversion to use
        if scan_date < datetime(20230208,'ConvertFrom','yyyymmdd')

            % Use original function

            % Check if mat file
            if strcmpi(p.subject_data.opm_noise{sub_num}(end-3:end), '.mat')

                % Convert the .mat data to BIDS format, suitable for
                % reading into FieldTrip
                convert_cMEG_mat_data(info);

            else

                % Convert the cMEG data to BIDS format, suitable for reading into
                % FieldTrip
                convert_cMEG_data(info);

            end

        else

            % Use new function
            convert_cMEG_data_tsv(info);

        end

        % Move files to their appropriate names
        movefile([p.directories.raw_dir(p.subject(sub_num), p.session(sub_num)) 'meg.bin'], p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), 'noise'));
        movefile([p.directories.raw_dir(p.subject(sub_num), p.session(sub_num)) 'meg.json'], p.raw_data.json(p.subject(sub_num), p.session(sub_num), p.run(sub_num), 'noise'));
        movefile([p.directories.raw_dir(p.subject(sub_num), p.session(sub_num)) 'positions.tsv'], p.raw_data.positions(p.subject(sub_num), p.session(sub_num), p.run(sub_num), 'noise'));
        movefile([p.directories.raw_dir(p.subject(sub_num), p.session(sub_num)) 'channels.tsv'], p.raw_data.channels(p.subject(sub_num), p.session(sub_num), p.run(sub_num), 'noise'));
        
        fprintf('done\n');
        
    end
    
    %% COPY THE SENS2MRI
    
    if ~exist(p.raw_data.sens2mri(p.subject(sub_num), p.session(sub_num)), 'file') || clobber == 1
        
        fprintf('Copying sens2mri transform...');
    
        % Get the sens2mri transform
        d = strsplit(p.subject_data.opm_data{sub_num}, '/');
        prefix = strsplit(d{end}, '_');
        fname = [strjoin(d, '/') '/' strjoin(prefix(1:2), '_') '_sens2mri_transform.txt'];
        
        % Copy
        copyfile(fname, p.raw_data.sens2mri(p.subject(sub_num), p.session(sub_num)));
        
        fprintf('done\n');
        
    end
    
end

%% TRIAL CONFIGURATION

% Only run if doesn't exist, or if overwriting
if ~exist(p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1
    
    fprintf('   Epoching...');
    
    % If clobbering, delete file
    if exist(p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file')
        delete(p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
    end
    
    % Check whether epoching by .json file or by trial duration
    cfg_epoch = [];
    if strcmp(p.epoching.type, 'json')
        
        % Json
        cfg_epoch.trialdef = p.epoching.json;
        cfg_epoch.trialfun = @opm_trialfun;
        
    else
        
        % Duration
        cfg_epoch.trialdef.triallength = p.epoching.duration;
        cfg_epoch.event = [];
        
    end
    
    % Do epoching
    cfg_epoch.dataset = p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num));
    cfg_epoch.headerformat = 'opm_fil_mjt';
    cfg_epoch.dataformat = 'opm_fil_mjt';
    cfg_epoch.eventformat = 'opm_fil_mjt';
    cfg_epoch.feedback = 'no';
    cfg_epoch = ft_definetrial(cfg_epoch);
    
    % Run the function to get the behavioural data
    if strcmp(p.epoching.type, 'json')
        [~, behav_table] = opm_trialfun(cfg_epoch);
        cfg_epoch.behav_table = behav_table;
    end
    
    % Save fieldtrip configuration
    parsave_struct(p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), cfg_epoch, 'cfg');
    
    fprintf('%d trials.\n', size(cfg_epoch.trl, 1));
    
    fprintf('done\n');
    
end


