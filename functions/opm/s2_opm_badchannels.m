function s2_opm_badchannels(p, sub_num, clobber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s2_opm_badchannels(p, sub_num, clobber)
%
% Script to iteratively perform mean field coorection, plot PSD and
% sensor-level OPM data, and identify bad channels. Channels that are
% outliers for more than 2/3 of the frequencies are suggested to be
% removed.
%
% Inputs:
%   p                   (structure) Study configuration structure
%   sub_num             (ind)       Subject number to process
%   clobber             (logical)   Whether or not to overwrite existing
%                                   files
%
% Outputs:
%   "bad_channels"      (file)      File containing the names of the bad
%                                   channels
%
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION

% Iteratively run MFC and bad channel identification and mean field
% correction until there is no additional bad channels
if (~exist(p.badchannels.bad_channels(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || ~exist(p.badchannels.bad_channels_psd(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file')) || clobber == 1
    
    % Get rid of any bad channels that already exist (either
    % overwriting, or failed run)
    if exist(p.badchannels.bad_channels_last(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file')
        delete(p.badchannels.bad_channels_last(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)))
    end
    if exist(p.badchannels.bad_channels(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file')
        delete(p.badchannels.bad_channels(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)))
    end
    
    close all force;
    
    %% LOAD IN THE DATA
    
    % Read in the header
    hdr = ft_read_header(p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'headerformat', 'opm_fil_mjt');
    
    % Load in the MEG data, and demean
    cfg_loadds = [];
    cfg_loadds.dataset = p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num));
    cfg_loadds.headerformat = 'opm_fil_mjt';
    cfg_loadds.dataformat = 'opm_fil_mjt';
    cfg_loadds.eventformat = 'opm_fil_mjt';
    cfg_loadds.channel = hdr.grad.label;
    cfg_loadds.demean = 'yes';
    [~, data] = evalc('ft_preprocessing(cfg_loadds)');
    
    % Make sure the grad structure contains the same labels (to fix issues
    % where sensors were down)
    [~, ind] = intersect(data.grad.label, data.label, 'stable');
    f = fieldnames(data.grad);
    grad_orig = data.grad.label;
    for ii = 1:length(f)
        if length(data.grad.(sprintf('%s', f{ii}))) == length(grad_orig)
            data.grad.(sprintf('%s', f{ii})) = data.grad.(sprintf('%s', f{ii}))(ind, :);
        end
    end
    
    % Store the data
    data_orig = data;
    
    % Don't mean field correct yet
    
    %% BAD CHANNELS
    
    bad_chs_new = NaN;
    while ~isempty(bad_chs_new)
        
        % If bad channels have already been identified and we're not to
        % overwrite, load, otherwise, ask user to define
        if ~exist(p.badchannels.bad_channels_last(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file')
            
            % Plot the sensor-level data (if not auto)
            if ~p.badchannels.auto
                cfg_ts = [];
                if isfield(p.preproc, 'bpfreq')
                    cfg_psd.foi = p.preproc.bpfreq;
                end
                cfg_psd.trial_length = 1;
                cfg_psd.average = 'median';
                sensor_ts(data, cfg_ts);
            end
            
            % Calculate and plot the PSD (if not auto)
            cfg_psd = [];
            if isfield(p.preproc, 'bpfreq')
                cfg_psd.foi = p.preproc.bpfreq;
            end
            cfg_psd.channel = true;
            cfg_psd.trial_length = 1;
            cfg_psd.average = 'median';
            if ~p.badchannels.auto
                cfg_psd.plot = true;
            else
                cfg_psd.plot = false;
            end
            [pow_raw, freq] = custom_psd(data, cfg_psd);
            
            % Take the log transform
            pow = log10(pow_raw);
            
            % Only keep the PSD within the FOI
            pow = pow(freq >= cfg_psd.foi(1) & freq <= cfg_psd.foi(2), :, :);
            
            % Get rid of outlier frequencies (peaks from power
            % spectrum)
            pow_o = median(pow, 3);
            pow_o = pow_o(~isoutlier(median(pow_o, 2)), :);
            
            % Find outliers >75% of the time
            outliers = sum(isoutlier(pow_o, 2), 1)'./size(pow_o, 1);
            remove = data.label(outliers > 2/3);
            
            % Check if auto is specified
            if p.badchannels.auto
                % Use the suggestions
                bad_chs = split(cell2mat({strjoin(remove, ', ')}), {', ', ','})';
            else
                % Ask the user to manually enter the bad channels
                opts = [];
                opts.WindowStyle = 'normal';
                bad_ch_in = select_bad_channels(remove, gcf);
                bad_chs = bad_ch_in; %split(cell2mat(bad_ch_in), {', ', ','})';
            end
            writecell(bad_chs, p.badchannels.bad_channels_last(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
            
        else
            
            % Load in the bad channels
            bad_chs = readcell(p.badchannels.bad_channels_last(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
            
        end
        
        close all force;
        
        % Get original
        data = data_orig;
        
        % Eliminate the bad channels
        cfg = [];
        cfg.channel = setdiff(data.label, bad_chs);
        [~, data] = evalc('ft_selectdata(cfg, data)');
        
        % Also apply to grad
        [~, ind] = intersect(data.grad.label, data.label, 'stable');
        grad_orig = data.grad.label;
        f = fieldnames(data.grad);
        for ii = 1:length(f)
            if length(data.grad.(sprintf('%s', f{ii}))) == length(grad_orig)
                data.grad.(sprintf('%s', f{ii})) = data.grad.(sprintf('%s', f{ii}))(ind, :);
            end
        end
        
        % Interference supression
        if strcmp(p.interference.method, 'HMC')
            
            % Mean field correct
            N = data.grad.chanori;
            M = eye(length(N)) - N*pinv(N);
            data.trial{1} = M*data.trial{1};
            
        elseif strcmp(p.interference.method, 'tSSS')
            
            % Check if we can perform tSSS
            if ((p.interference.tsss.lin+1)^2+(p.interference.tsss.lout+1)^2-2) > size(data.grad.chanpos,1)
                error('Not enough channels remaining for tSSS!');
            end
            
            % Temporal signal space separation
            sss_outs = tSSS_opm(data.trial{1}, data.grad.chanpos, data.grad.chanori, p.interference.tsss.lin, p.interference.tsss.lout, 1, p.interference.tsss.corr);
            data.trial{1} = sss_outs.tSSS;
            
        else
            
            % Error out
            error('Specify either HMC or tSSS for the interference supression method!');
            
        end
        
        % Plot the sensor-level data
        if ~p.badchannels.auto
            cfg_ts = [];
            if isfield(p.preproc, 'bpfreq')
                cfg_ts.foi = p.preproc.bpfreq;
            end
            cfg_ts.trial_length = 1;
            cfg_ts.average = 'median';
            sensor_ts(data, cfg_ts);
        end
        
        % Plot the PSD
        cfg_psd = [];
        if isfield(p.preproc, 'bpfreq')
            cfg_psd.foi = p.preproc.bpfreq;
        end
        cfg_psd.channel = true;
        cfg_psd.trial_length = 1;
        cfg_psd.average = 'median';
        if ~p.badchannels.auto
            cfg_psd.plot = true;
        else
            cfg_psd.plot = false;
        end
        [pow, freq] = custom_psd(data, cfg_psd);
        
        % Only keep the PSD within the FOI
        pow = pow(freq >= cfg_psd.foi(1) & freq <= cfg_psd.foi(2), :, :);
        
        % Get rid of outlier frequencies (peaks from power
        % spectrum)
        pow_o = median(pow, 3);
        pow_o = pow_o(~isoutlier(median(pow_o, 2)), :);
        
        % Find outliers >75% of the time
        outliers = sum(isoutlier(pow_o, 2), 1)'./size(pow_o, 1);
        remove = data.label(outliers > 2/3);
        
        % Check if auto is specified
        if p.badchannels.auto
            % Use the suggestions
            bad_chs_new = split(cell2mat({strjoin(remove, ', ')}), {', ', ','})';
        else
            % Check if any more bad channels need to be removed
            opts = [];
            opts.WindowStyle = 'normal';
            bad_ch_in = select_bad_channels(remove, gcf);
            bad_chs_new = bad_ch_in; %split(cell2mat(bad_ch_in), {', ', ','})';
        end
        
        % Check if any new
        if ~isempty(cell2mat(bad_chs_new))
            
            % Concatenate with existing bad channels
            if ~cellfun(@isempty, bad_chs)
                bad_chs = cat(2, bad_chs, bad_chs_new);
            else
                bad_chs = bad_chs_new;
            end
            writecell(bad_chs, p.badchannels.bad_channels_last(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
            
        else
            
            % Exit loop
            bad_chs_new = [];
            
            % Save bad channels
            copyfile(p.badchannels.bad_channels_last(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), p.badchannels.bad_channels(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
            delete(p.badchannels.bad_channels_last(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
            
        end
        
        
    end
    
    %% CHECK NUMBER OF BAD CHANNELS
    
    % Check if we can perform tSSS
    if strcmp(p.interference.method, 'tSSS')
        if ((p.interference.tsss.lin+1)^2+(p.interference.tsss.lout+1)^2-2) > length(setdiff(data_orig.label, bad_chs))
            error('Not enough channels remaining for tSSS!');
        end
    end
    
    %% MAKE FINAL PSD FIGURE
    
    close all force;
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    fig.Name = 'PSD';
    t = tiledlayout(fig, 1, 2);
    
    % Plot the original
    nexttile
    
    % Get original
    data = data_orig;
    
    % Interference supression
    if strcmp(p.interference.method, 'HMC')
        % Mean field correct
        N = data.grad.chanori;
        M = eye(length(N)) - N*pinv(N);
        data.trial{1} = M*data.trial{1};
    elseif strcmp(p.interference.method, 'tSSS')
        % Temporal signal space separation
        sss_outs = tSSS_opm(data.trial{1}, data.grad.chanpos, data.grad.chanori, p.interference.tsss.lin, p.interference.tsss.lout, 1, p.interference.tsss.corr);
        data.trial{1} = sss_outs.tSSS;
    else
        % Error out
        error('Specify either HMC or tSSS for the interference supression method!');
    end
    
    % Plot the PSD
    cfg_psd = [];
    if isfield(p.preproc, 'bpfreq')
        cfg_psd.foi = p.preproc.bpfreq;
    end
    cfg_psd.channel = true;
    cfg_psd.trial_length = 1;
    cfg_psd.average = 'median';
    cfg_psd.plot = true;
    custom_psd(data, cfg_psd, fig);
    
    % Plot the new
    nexttile
    
    % Get original
    data = data_orig;
    
    % Eliminate the bad channels
    cfg = [];
    cfg.channel = setdiff(data.label, bad_chs);
    [~, data] = evalc('ft_selectdata(cfg, data)');
    
    % Also apply to grad
    [~, ind] = intersect(data.grad.label, data.label, 'stable');
    grad_orig = data.grad.label;
    f = fieldnames(data.grad);
    for ii = 1:length(f)
        if length(data.grad.(sprintf('%s', f{ii}))) == length(grad_orig)
            data.grad.(sprintf('%s', f{ii})) = data.grad.(sprintf('%s', f{ii}))(ind, :);
        end
    end
    
    % Interference supression
    if strcmp(p.interference.method, 'HMC')
        % Mean field correct
        N = data.grad.chanori;
        M = eye(length(N)) - N*pinv(N);
        data.trial{1} = M*data.trial{1};
    elseif strcmp(p.interference.method, 'tSSS')
        % Temporal signal space separation
        sss_outs = tSSS_opm(data.trial{1}, data.grad.chanpos, data.grad.chanori, p.interference.tsss.lin, p.interference.tsss.lout, 1, p.interference.tsss.corr);
        data.trial{1} = sss_outs.tSSS;
    else
        % Error out
        error('Specify either HMC or tSSS for the interference supression method!');
    end
    
    % Plot the PSD
    cfg_psd = [];
    if isfield(p.preproc, 'bpfreq')
        cfg_psd.foi = p.preproc.bpfreq;
    end
    cfg_psd.channel = true;
    cfg_psd.trial_length = 1;
    cfg_psd.average = 'median';
    cfg_psd.plot = true;
    custom_psd(data, cfg_psd, fig);
    
    % Set the ylimits to be the same
    yl = [min(t.Children(1).YLim(1), t.Children(2).YLim(1)), max(t.Children(1).YLim(2), t.Children(2).YLim(2))];
    t.Children(1).YLim = yl;
    t.Children(2).YLim = yl;
    
    % Set the titles
    t.Children(2).Title.String = 'Original';
    t.Children(1).Title.String = p.interference.method;
    
    % Save
    saveas(gcf, p.badchannels.bad_channels_psd(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
    close all
    
    
end

function bad = select_bad_channels(preselected, fig)
% SELECT_BAD_CHANNELS  Click traces in an existing figure to toggle bad channels;
% press Enter to finish. Channel names are read from each line's datatip.
%   bad = select_bad_channels(preselected)                 % uses gcf
%   bad = select_bad_channels(preselected, fig)
%   bad = select_bad_channels(preselected, fig, rowName)   % e.g. 'Channel', if the
%                                                          % datatip has several text rows
%   Returns a cellstr of bad channel names.

if nargin < 1, preselected = {}; end
if nargin < 2 || isempty(fig), fig = gcf; end
preselected = cellstr(preselected);
labelRowNumber = 3;

% ---------- read channel names from the datatips ----------
allLns  = findobj(fig, 'Type', 'line');
lineLab = cell(numel(allLns), 1);
for lineIdx = 1:numel(allLns)
    rows = allLns(lineIdx).DataTipTemplate.DataTipRows;
    if numel(rows) >= labelRowNumber
        lineLab{lineIdx} = rows(labelRowNumber).Label;
    end
end
keep    = ~cellfun(@isempty, lineLab);
lns     = allLns(keep);
lineLab = lineLab(keep);
if isempty(lns)
    error('No lines with a text datatip row were found in this figure.');
end
labels     = unique(lineLab, 'stable')';
[~, chIdx] = ismember(lineLab, labels);      % line k belongs to channel chIdx(k)

isBad       = ismember(labels, preselected);
extra       = setdiff(preselected, labels, 'stable');   % preselected but not in figure
lastClicked = '';
figClosed   = false;

% ---------- remember original state so it can be restored ----------
origColor = arrayfun(@(l) l.Color, lns, 'UniformOutput', false);
origWidth = [lns.LineWidth];
oldKey    = fig.KeyPressFcn;
oldClose  = fig.CloseRequestFcn;
axList    = findobj(fig, 'Type', 'axes');

for a = axList', disableDefaultInteractivity(a); end
set(lns, 'ButtonDownFcn', @onClick);
fig.KeyPressFcn     = @onKey;
fig.CloseRequestFcn = @onCloseFig;

banner = annotation(fig, 'textbox', [0 0.93 1 0.07], 'EdgeColor', 'none', ...
                    'BackgroundColor', [1 1 0.85], 'FaceAlpha', 0.9, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                    'Interpreter', 'none', 'FontSize', 9);

redraw();
figure(fig);                                 % bring to front so it receives the Enter key
uiwait(fig);

bad = [labels(isBad), extra];

% ---------- clean up ----------
if isvalid(fig)
    if figClosed
        delete(fig);
    else
        delete(banner);
        fig.KeyPressFcn     = oldKey;
        fig.CloseRequestFcn = oldClose;
        set(lns(isvalid(lns)), 'ButtonDownFcn', '');
        for a = axList', enableDefaultInteractivity(a); end
    end
end

% ---------- callbacks ----------
    function onClick(src, ~)
        c = chIdx(find(lns == src, 1));
        isBad(c) = ~isBad(c);
        lastClicked = labels{c};
        redraw();
    end

    function onKey(~, evt)
        switch evt.Key
            case {'return', 'enter'}
                uiresume(fig);
            case {'c', 'delete', 'backspace'}
                isBad(:) = false;   % unmark every channel in the figure
                extra = {};         % also drop preselected channels not in the figure
                lastClicked = '';
                redraw();
        end
    end

    function onCloseFig(~, ~)
        figClosed = true;
        uiresume(fig);
    end

% ---------- helper ----------
    function redraw()
        for k = 1:numel(lns)
            if isBad(chIdx(k))
                set(lns(k), 'Color', [0.85 0.1 0.1], 'LineWidth', 2);
                uistack(lns(k), 'top');
            else
                set(lns(k), 'Color', origColor{k}, 'LineWidth', origWidth(k));
            end
        end
        banner.String = {sprintf(['Click traces to toggle bad channels, C to clear all, ' ...
            'Enter when done.   Last clicked: %s'], lastClicked), ...
            ['Bad: ' strjoin([labels(isBad), extra], ', ')]};
    end
end

function lab = tipLabel(l, rowName)
% Return the first text value from a line's DataTipTemplate ('' if none).
lab = '';
try
    rows = l.DataTipTemplate.DataTipRows;
catch
    return
end
for r = rows(:)'
    if ~isempty(rowName) && ~strcmpi(r.Label, rowName), continue; end
    v = r.Value;
    if iscell(v) || isstring(v) || iscategorical(v)
        lab = char(string(v(1)));
        return
    end
end
end


end