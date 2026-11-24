function convert_cMEG_data(cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_cMEG_data(info)
%
% Read OPM data from Cerca in a BIDS structure suitable for conversion to
% FieldTrip. This function applys the conversion factor recorded by the 
% operator (Gain x1 has conversion 2.7 V/nT), and is applicable to v1 cMEG
% data
%
% Inputs:
%   info                (structure) Input structure
%   info.filename       (string)    Path of cerca data file
%   info.pos_dir        (string)    Path containing helmet layouts
%   info.output_dir     (string)    Path of directory to save output
%
% Outputs:
%   "output_dir"        (directory)
%       "channels.tsv"  (file)      BIDS channels table
%       "meg.bin"       (file)      FieldTrip OPM data
%       "meg.json"      (file)      BIDS sidecar JSON
%       "positions.tsv" (file)      FieldTrip OPM sensor positions
%
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHECK INPUTS

% Make sure all the required inputs exist
if ~isfield(cfg, 'filename')
    error('Must define filename in the configuration structure');
end
if ~isfield(cfg, 'pos_dir')
    error('Must define pos_dir in the configuration structure');
end
if ~isfield(cfg, 'output_dir')
    error('Must define output_dir in the configuration structure');
end

%% READ IN RAW DATA

% Extract filename
if strcmpi(num2str(cfg.filename(end-4:end)),'.cMEG')
    filename = cfg.filename(1:end-5);
else
    folders = dir([cfg.filename '/*.cMEG']);
    if length(folders) > 1
        error('Multiple .cMEG files found! Please specify the exact .cMEG file instead');
    else
        filename = [folders.folder '/' folders.name(1:end-5)];
    end
end

% Open the file
fname = [filename '.cMEG'];
fid = fopen(fname,'r','ieee-be');
finfo = dir(fname);
fsize = finfo.bytes;
Adim_conv = [2^32; 2^16; 2^8; 1]; % Array dimension conversion table

% Preallocation
I = 0;
while fsize ~= ftell(fid)
    dims = [];
    for n = 1:2 % 2 dimensions (Nch x Time)
        temp = fread(fid,4);
        temp = sum(Adim_conv.*temp);
        dims = cat(1, dims, temp);
    end
    I = I + 1;
    temp = fread(fid,prod(dims),'double',0,'ieee-be');  % Skip the actual data (much faster for some reason)
end
fseek(fid,0,'bof'); % Reset the cursor
data1 = repmat({NaN*ones(dims')},I,1);  % Preallocate space, assumes each section is the same
zeros(size(temp)); % Just to get rid of warning

% Load and parse data
for j = 1:I
    dims = [];  % Actual array dimensions
    for n = 1:2 % 2 dimensions (Nch x Time)
        temp = fread(fid,4);
        temp = sum(Adim_conv.*temp);
        dims = cat(1, dims, temp);
    end
    clear temp n
    temp = fread(fid,prod(dims),'double',0,'ieee-be');
    
    % Reshape the data into the correct array configuration
    for n = 1:2
        temp = reshape(temp,dims(2-n+1),[])';
    end
    data1{j} = temp;  % Save the data
end
fclose(fid);  % Close the file

% Reshape into sensible order
data_orig = NaN*zeros(size(data1{1},2),size(data1{1},1).*size(data1,1));
for n = 1:size(data1,1)
    %     disp(100*n/length(data1))
    data_range = (n-1)*size(data1{1},1)+1:n*size(data1{1},1);
    data_orig(:,data_range) = data1{n}'; % data in Nch+triggers+1 x time, but channel one is the time
end

% Get sampling frequency
samp_frequency = round(1/(data_orig(1,2)-data_orig(1,1)));

% Get rid of first sensor (noise?)
data_orig = data_orig(2:end, :);

%% READ IN SESSION INFO

% Read session info
Session_info_file = [filename '_SessionInfo.txt'];
fidSI = fopen(Session_info_file);
finfoSI = dir(Session_info_file);
fsizeSI = finfoSI.bytes;
count = 0;
while fsizeSI ~= ftell(fidSI)
    count = count + 1;
    Session_info{count,1} = fgetl(fidSI);
end
fclose(fidSI);

% Get channel names
for n = 1:size(Session_info,1)
    if strfind(Session_info{n},'Sensor Names:')
        Chan_info_string_id = n;
    end
end
Chan_info_string = Session_info{Chan_info_string_id};
cn1 = strsplit(Chan_info_string,', ');
cn11 = strsplit(cn1{1},':  ');
chan_names = [cn11(2) cn1(2:end)]';

%% SEPARATE TRIGGER AND OPM DATA

% Extract the triggers
tIdx = strfind(chan_names,'Trigger');
tIdx = find(~cellfun(@isempty, tIdx));
triggers = data_orig(tIdx,:);

% Extract the OPM data
opm_data_orig = data_orig;
opm_data_orig(tIdx,:) = [];

% Get the name of the OPM channels and the triggers
opm_names = chan_names;
trigger_names = opm_names(tIdx,:);
opm_names(tIdx,:) = [];

% Split the name from the axis
names_only = cell(size(opm_names, 1), 1);
axis_only = cell(size(opm_names, 1), 1);
for n = 1:size(opm_names,1)
    tmp = strsplit(opm_names{n});
    names_only(n,1) = tmp(1);
    axis_only(n,1) = tmp(2);
end 

%% ADJUST OPM DATA

% Extract the OPM gain
gain_line = split(Session_info{contains(Session_info, 'OPM Gain:')},': ');
gain_conv = 2.7.*str2double(gain_line{2}(1:end-1));

% Apply conversion factor and convert to T
opm_data_orig = (1e-9.*opm_data_orig)./(gain_conv);

% Load the sensor order
sensor_order = load([fname(1:end-5) '_sensor_order.mat']);

% Iterate through reorder the OPM data in layout order according to the
% sensor order
loc_info = cell(size(sensor_order.T.Name));
opm_data = zeros(size(opm_data_orig));
sensor_names = cell(size(opm_data_orig, 1), 1);
opm_data_mat = [];
count = 0;
for ii = 1:size(sensor_order.T, 1)
    if ~isempty(sensor_order.T.Name{ii})
        sens_name = strsplit(sensor_order.T.Name{ii});
        trigger_sensors = find(strcmpi(sens_name(1), names_only));
        for m = 1:length(trigger_sensors)
            count = count + 1;
            opm_data_mat.(['Sensor_' sens_name{1}]).(axis_only{trigger_sensors(m)}(2)) = opm_data_orig(trigger_sensors(m),:);
            loc_info{ii} = ['Sensor ' sens_name{1}];
            sensor_names{count} = [sens_name{1} ' ' axis_only{trigger_sensors(m)}(2)];
            opm_data(count, :) = opm_data_orig(trigger_sensors(m),:);
        end
    end
end

% Load in the helmet info here
helmet = split(Session_info(contains(Session_info, 'Helmet:') & ~contains(Session_info, 'too Large')), ': ');
helmet = strrep(helmet{2}, ' ', '_');
Helmet_info = load([cfg.pos_dir '/' helmet '.mat']);
Helmet_info = Helmet_info.Helmet_info;
sens_pos = zeros(size(sensor_names, 1), 3);
sens_ors = zeros(size(sensor_names, 1), 3);
for ii = 1:size(sensor_names,1)
    if ~isempty(sensor_names{ii})
        tmp = split(sensor_names{ii});
        loc_idx = unique(find(strcmpi(['Sensor ' tmp{1}], loc_info)));
        sens_pos(ii,:) = Helmet_info.sens_pos(loc_idx,:);
        sens_ors(ii,:) = Helmet_info.(['sens_ors_' tmp{2}])(loc_idx,:);
    end
end

%% SAVE BIDS FILES

% Make sure BIDS directory exists
if ~exist(cfg.output_dir, 'dir')
    mkdir(cfg.output_dir);
end

% Create the channels structure and write
channels = [];
channels.name = [sensor_names; trigger_names];
channels.type = [repmat({'MEGMAG'},size(sensor_names,1), 1); repmat({'TRIG'}, size(triggers, 1), 1)];
channels.units = [repmat({'fT'}, size(sensor_names,1), 1); repmat({'V'}, size(trigger_names, 1), 1)];
channels.status = repmat({'good'}, size(channels.name,1), 1);
writetable(struct2table(channels), [cfg.output_dir '/channels.tsv'], 'Delimiter', '\t', 'FileType', 'text');

% Write data file
fid = fopen([cfg.output_dir '/meg.bin'], 'w');
fwrite(fid, [opm_data; triggers], 'double', 'b');
fclose(fid);

% Write header file
meg=[];
meg.SamplingFrequency = samp_frequency;
meg.Helmet = helmet;
meg = jsonencode(meg);
fid = fopen([cfg.output_dir '/meg.json'], 'w');
fprintf(fid, meg);
fclose(fid);

% Create the positions structure and write
positions=[];
positions.name= sensor_names;
positions.Px = sens_pos(:,1);
positions.Py = sens_pos(:,2);
positions.Pz = sens_pos(:,3);
positions.Ox = sens_ors(:,1);
positions.Oy = sens_ors(:,2);
positions.Oz = sens_ors(:,3);
writetable(struct2table(positions), [cfg.output_dir '/positions.tsv'], 'Delimiter', '\t', 'FileType', 'text');

end