function cMEG_combine_v28(all_files, savename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cMEG_combine_v28(all_files, savename)
%
% Combined OPM cMEG files in the order of the "all_files" variable and
% stores the combined output in "savename"
%
% Inputs:
%   all_files	(cell)      Cell array of strings specifying file locations
%   savename    (string)    Save name
% 
% Outputs:
%   "savename"  (file)      cMEG file containing combined data
%
% Created by: Marlee M Vandewouw
% Created on: March 29th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Iterate through the files
data_all = [];
for num_fname = 1:length(all_files)
    
    % Clear all variables from last iteration except these ones
    clearvars -except filename_base all_files num_fname data_all savename
    
    % Filename of current file
    fname = [all_files{num_fname}];
    
    % Open the file
    fid = fopen(fname,'rb','ieee-be');
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
        dims = []; % Actual array dimensions
        for n = 1:2 % 2 dimensions (Nch x Time)
            temp = fread(fid,4);
            temp = sum(Adim_conv.*temp);
            dims = cat(1, dims, temp);
        end
        clear temp n
        temp = fread(fid,prod(dims),'double',0,'ieee-be');
        data1{j} = temp;  % Save the data
    end
    fclose(fid);  % Close the file
    
    % Concatenate
    data_all = cat(1, data_all, data1);
end

% Extract the number of chanels and number of sample segs
Nch = dims(1);
Nsamp_seg = dims(2);

% Write combined data
fid = fopen(savename,'w','ieee-be');
for n_data = 1:size(data_all,1)
    rder = Nch;
    rder2 = Nsamp_seg;
    Nch_vec = zeros(length(Adim_conv), 1);
    Nsamp_vec = zeros(length(Adim_conv), 1);
    for nn = 1:length(Adim_conv)
        val = floor(rder./Adim_conv(nn));
        rder = rder - Adim_conv(nn).*val;
        Nch_vec(nn) = val;
        val2 = floor(rder2./Adim_conv(nn));
        rder2 = rder2 - Adim_conv(nn).*val2;
        Nsamp_vec(nn) = val2;
    end
    hdr = [Nch_vec, Nsamp_vec];
    % Write header
    fwrite(fid,hdr);
    % Write data segment
    fwrite(fid,data_all{n_data},'double');
end
fclose(fid);

end