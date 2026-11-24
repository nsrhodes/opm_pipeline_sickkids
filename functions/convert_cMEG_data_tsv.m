function convert_cMEG_data_tsv(info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_cMEG_data_v2p28(info)
%
% Read OPM data from Cerca in a BIDS structure suitable for conversion to
% FieldTrip. This function applys the conversion factor recorded by the 
% operator, and is applicable to v2 cMEG data
%
% Inputs:
%   info                    (structure) Input structure
%   info.filename           (string)    Path of cerca data file
%   info.pos_dir            (string)    Path containing helmet layouts
%   info.output_dir         (string)    Path of directory to save output
%
% Outputs:
%   "output_dir"            (directory)
%       "channels.tsv"      (file)      BIDS channels table
%       "meg.bin"           (file)      FieldTrip OPM data
%       "meg.json"          (file)      BIDS sidecar JSON
%       "positions.tsv"     (file)      FieldTrip OPM sensor positions
%       "Layout_Info.mat"   (file)      MATLAB structure containing
%                                       information on the sensor layout
%
% Created by: Marlee M Vandewouw
% Created on: March 29th,, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHECK INPUTS

% Make sure all the required inputs exist
if ~isfield(info, 'filename')
    error('Must define filename in the configuration structure');
end
if ~isfield(info, 'pos_dir')
    error('Must define pos_dir in the configuration structure');
end
if ~isfield(info, 'output_dir')
    error('Must define output_dir in the configuration structure');
end

%% GET ALL CMEG FILES AND COMBINE

% Extract filename
if strcmpi(info.filename(end-4:end),'.cMEG')
    splitted_filename = split(info.filename, '/');
    foldername = char(join(splitted_filename(1:end-1), '/'));
else
    foldername = info.filename;
end

% Get all the subdirectories
cMEG_files = dir([foldername '/*.cMEG']);

% Check if combined file has already been made
ind = find(contains({cMEG_files.name}, '_meg.cMEG'));

% If the combined file hasn't been made - make a temporary one in the
% output directory
if isempty(ind)
    
    % Create cell array of all the files to combine
    all_files = cell(length(cMEG_files), 1);
    for ii = 1:length(all_files)
        all_files{ii} = [foldername '/' cMEG_files(ii).name];
    end
    
    % Combine
    cMEG_file = [info.output_dir '/temp_combined.cMEG'];
    cMEG_combine_v28(all_files, cMEG_file)
        
else

    % Use the combined file
    cMEG_file = info.filename;%[foldername '/' cMEG_files(ind).name];
    
end

%% READ IN RAW DATA

% Open the file
fid = fopen(cMEG_file,'r','ieee-be');
finfo = dir(cMEG_file);
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
    dims = []; % Actual array dimensions
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
    data_range = (n-1)*size(data1{1},1)+1:n*size(data1{1},1);
    data_orig(:,data_range) = data1{n}'; % data in Nch+triggers+1 x time, but channel one is the time
end

% Get the data in volts
data_orig = data_orig(2:end, :);

%% READ IN CHANNEL INFORMATION

% Find the channels file
channel_info_file = dir([foldername '/*_channels.tsv']);
channel_info_file = [channel_info_file.folder '/' channel_info_file.name];
channel_info  = tdfread(channel_info_file);

% Fix the strings
CI_fields = fields(channel_info);
for n = 1:size(CI_fields, 1)
    try
        channel_info.(CI_fields{n}) = cellstr(channel_info.(CI_fields{n}));
    end
end

%% READ IN SESSION INFO

% Read session info
sesh_info_file = dir([foldername '/*_meg.json']);
sesh_fid = fopen([sesh_info_file.folder '/' sesh_info_file.name]);
raw = fread(sesh_fid,inf);
str = char(raw');
fclose(sesh_fid);
Session_info = jsondecode(str);

%% SEPARATE TRIGGER AND OPM DATA

% Extract the triggers
tIdx = strfind(channel_info.name, 'Trigger');
tIdx = find(~cellfun(@isempty, tIdx));
triggers = data_orig(tIdx,:);
trigger_names = channel_info.name(tIdx);

% Extract the OPM data
opm_data_orig = data_orig;
opm_data_orig(tIdx,:) = [];

% Get the new channel info
oIdx = 1:length(channel_info.name);
oIdx(tIdx) = [];
f = fieldnames(channel_info);
for ii = 1:length(f)
    channel_info.(sprintf('%s', f{ii})) = channel_info.(sprintf('%s', f{ii}))(oIdx, :);
end

%% ADJUST OPM DATA

% Iterate over the sensors
data = zeros(size(opm_data_orig));
for nn = 1:size(data, 1)
    
    % Extract the OPM gain
    gain_conv = str2double(channel_info.nT0x2FV(nn, :));
    if isnan(gain_conv)
        gain_conv = 1;
    end

    % Apply conversion factor and convert to T
    data(nn, :) = (1e-9.*opm_data_orig(nn, :))./(gain_conv);
    
end

% Load in the layout 
layout_info_file = dir([foldername '/*_HelmConfig.tsv']);
Layout_Info = tdfread([layout_info_file.folder '/' layout_info_file.name]);
last_line_idx = strmatch('Helmet:',Layout_Info.Name);
field_names = fieldnames(Layout_Info);
for n = 1:size(field_names,1)
    Layout_Info.(field_names{n})(last_line_idx,:) = [];
    if max(strcmpi(field_names{n},{'Px';'Py';'Pz';'Ox';'Oy';'Oz';'Layx';'Layz'}))
        if ischar(Layout_Info.(field_names{n})(1))
            Layout_Info.(field_names{n}) = str2double(cellstr(Layout_Info.(field_names{n})));
        end
    end
end

% Fix the names
Layout_Info.Name = cellstr(Layout_Info.Name);

% Get the slot/sensor pairs
Layout_Info.SlotSensorPairs = [cellstr(Layout_Info.Name) cellstr(Layout_Info.Sensor)];
CI_name = channel_info.name;
for n = 1:size(CI_name,1)
    CI_name{n} = CI_name{n}(isstrprop(CI_name{n},'alphanum'));
end
for n = 1:size(Layout_Info.SlotSensorPairs,1)
    sens_name = Layout_Info.SlotSensorPairs{n,2};
    sens_name = sens_name(isstrprop(sens_name,'alphanum'));

    data_idx = find(strcmpi(sens_name,CI_name));
    if isempty(data_idx)
        Layout_Info.SlotSensorPairs{n,3} = 0;
    else
        Layout_Info.SlotSensorPairs{n,3} = data_idx;
    end
end

% Get the sensor positions and orientations
sens_pos = zeros(length(channel_info.name), 3);
sens_ors = zeros(length(channel_info.name), 3);
for ii = 1:length(channel_info.name)
    ind = find(strcmp(Layout_Info.SlotSensorPairs(:, 2), strrep(strrep(channel_info.name{ii}, '[', ''), ']', '')));
    if ~isempty(ind)
        sens_pos(ii, :) = [Layout_Info.Px(ind), Layout_Info.Py(ind), Layout_Info.Pz(ind)];
        sens_ors(ii, :) = [Layout_Info.Ox(ind), Layout_Info.Oy(ind), Layout_Info.Oz(ind)];
    end
end

%% SAVE BIDS FILES

% Create the channels structure and write
sensor_names = channel_info.name;
channels = [];
channels.name = [sensor_names; trigger_names];
channels.type = [repmat({'MEGMAG'},size(sensor_names,1), 1); repmat({'TRIG'}, size(triggers, 1), 1)];
channels.units = [repmat({'fT'}, size(sensor_names,1), 1); repmat({'V'}, size(trigger_names, 1), 1)];
channels.status = repmat({'good'}, size(channels.name,1), 1);
writetable(struct2table(channels), [info.output_dir '/channels.tsv'], 'Delimiter', '\t', 'FileType', 'text');

% Write data file
fid = fopen([info.output_dir '/meg.bin'], 'w');
fwrite(fid, [data; triggers], 'double', 'b');
fclose(fid);

% Write header file
meg = jsonencode(Session_info);
fid = fopen([info.output_dir '/meg.json'], 'w');
fprintf(fid, meg);
fclose(fid);

% Create the positions structure and write
positions=[];
positions.name = channel_info.name;
positions.Px = sens_pos(:,1);
positions.Py = sens_pos(:,2);
positions.Pz = sens_pos(:,3);
positions.Ox = sens_ors(:,1);
positions.Oy = sens_ors(:,2);
positions.Oz = sens_ors(:,3);
writetable(struct2table(positions), [info.output_dir '/positions.tsv'], 'Delimiter', '\t', 'FileType', 'text');

% Save the layout info
parsave_struct([info.output_dir '/Layout_Info.mat'], Layout_Info, 'Layout_Info'); 

% Delete temporary combined file
if isempty(ind)
    delete(cMEG_file);
end

end