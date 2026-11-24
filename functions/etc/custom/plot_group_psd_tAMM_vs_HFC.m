function plot_group_psd_tAMM_vs_HFC(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_group_psd_tAMM_vs_HFC(p)
%
% Plots group-level PSDs comparing tAMM and HFC interference suppression,
% after bad channels have been identified via s2_opm_badchannels.
%
% One line per participant is drawn on each panel (channel-averaged PSD),
% with both panels sharing the same y-axis limits.
%
% Inputs:
%   p   (structure)  Study configuration structure
%
% Outputs:
%   Figure saved to p.badchannels.group_psd_comparison (define in p)
%
% Requires: s2_opm_badchannels already run for all subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_subs = size(p.subject_data, 1);

% Pre-allocate cell arrays to hold per-subject PSDs and freq axis
psd_tAMM = cell(n_subs, 1);
psd_HFC  = cell(n_subs, 1);
freq_all = cell(n_subs, 1);

%% LOOP OVER SUBJECTS

for ss = 1:n_subs

    fprintf('Processing subject: %s, session: %s, task: %s, run: %s\n', ...
        p.subject(ss), p.session(ss), p.task(ss), p.run(ss));

    % ------------------------------------------------------------------ %
    %  LOAD RAW DATA
    % ------------------------------------------------------------------ %

    hdr = ft_read_header( ...
        p.raw_data.bin(p.subject(ss), p.session(ss), p.run(ss), p.task(ss)), ...
        'headerformat', 'opm_fil_mjt');

    cfg_loadds                = [];
    cfg_loadds.dataset        = p.raw_data.bin(p.subject(ss), p.session(ss), p.run(ss), p.task(ss));
    cfg_loadds.headerformat   = 'opm_fil_mjt';
    cfg_loadds.dataformat     = 'opm_fil_mjt';
    cfg_loadds.eventformat    = 'opm_fil_mjt';
    cfg_loadds.channel        = hdr.grad.label;
    cfg_loadds.demean         = 'yes';
    [~, data] = evalc('ft_preprocessing(cfg_loadds)');

    % Fix grad to match loaded channels
    [~, ind] = intersect(data.grad.label, data.label, 'stable');
    data = fix_grad(data, ind);

    % ------------------------------------------------------------------ %
    %  REMOVE BAD CHANNELS
    % ------------------------------------------------------------------ %

    bad_ch_file = p.badchannels.bad_channels( ...
    p.subject(ss), p.session(ss), p.run(ss), p.task(ss));
    
    if ~exist(bad_ch_file, 'file')
        warning('Bad channel file not found for subject %s — skipping.', p.subject(ss));
        continue
    end
    
    % Read raw text and split on commas/newlines only (preserves channel
    % names that contain spaces, e.g. '99 [Y]')
    raw     = strsplit(strtrim(fileread(bad_ch_file)), {',', '\n', '\r'});
    bad_chs = strtrim(raw);
    bad_chs = bad_chs(~cellfun(@isempty, bad_chs));

    cfg_sel         = [];
    cfg_sel.channel = setdiff(data.label, bad_chs);
    [~, data_clean] = evalc('ft_selectdata(cfg_sel, data)');

    % Fix grad after channel removal
    [~, ind] = intersect(data_clean.grad.label, data_clean.label, 'stable');
    data_clean = fix_grad(data_clean, ind);

    % ------------------------------------------------------------------ %
    %  PSD CONFIG (shared)
    % ------------------------------------------------------------------ %

    cfg_psd              = [];
    cfg_psd.trial_length = 1;
    cfg_psd.average      = 'median';
    cfg_psd.channel      = true;
    cfg_psd.plot         = false;
    if isfield(p.preproc, 'bpfreq')
        cfg_psd.foi = p.preproc.bpfreq;
    end

    % ------------------------------------------------------------------ %
    %  tAMM
    % ------------------------------------------------------------------ %

    data_tAMM = data_clean;  % work on a copy
    amm_outs  = tAMM_opm( ...
        data_tAMM.trial{1}, ...
        data_tAMM.grad.chanpos, ...
        data_tAMM.grad.chanori, ...
        p.interference.tamm.lin, ...
        p.interference.tamm.lout, ...
        p.interference.tamm.corr);
    data_tAMM.trial{1} = amm_outs.tAMM_data;

    [pow_tAMM, freq] = custom_psd(data_tAMM, cfg_psd);

    % ------------------------------------------------------------------ %
    %  HFC  (homogeneous field correction / mean field correction)
    % ------------------------------------------------------------------ %

    data_HFC = data_clean;  % fresh copy — tAMM result discarded
    N = data_HFC.grad.chanori;
    M = eye(length(N)) - N * pinv(N);
    data_HFC.trial{1} = M * data_HFC.trial{1};

    [pow_HFC, ~] = custom_psd(data_HFC, cfg_psd);

    % ------------------------------------------------------------------ %
    %  AVERAGE ACROSS CHANNELS, LOG-TRANSFORM, TRIM TO FOI
    % ------------------------------------------------------------------ %

    foi_mask = freq >= cfg_psd.foi(1) & freq <= cfg_psd.foi(2);
    freq_foi = freq(foi_mask);

    % pow is [freq x channels x trials] — median across trials already
    % applied inside custom_psd; average across channels here
    % Trim to FOI
    pow_foi_tAMM = pow_tAMM(foi_mask, :, :);   % [freq x channels x epochs]
    pow_foi_HFC  = pow_HFC(foi_mask,  :, :);
    
    % Median across epochs, then mean across channels -> [freq x 1]
    psd_tAMM{ss} = mean(log10(median(pow_foi_tAMM, 3)), 2);
    psd_HFC{ss}  = mean(log10(median(pow_foi_HFC,  3)), 2);
    freq_all{ss} = freq(foi_mask);

end

%% PLOT

fig = figure('units', 'normalized', 'outerposition', [0 0 1 0.6]);
fig.Name = 'Group PSD: tAMM vs HFC';
t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Colour map — one colour per subject
cmap = lines(n_subs);

ax_tAMM = nexttile; hold on;
ax_HFC  = nexttile; hold on;

legend_entries = {};

for ss = 1:n_subs

    if isempty(psd_tAMM{ss}), continue; end   % subject was skipped

    col  = cmap(ss, :);
    subj_label = sprintf('%s', p.subject(ss));

    plot(ax_tAMM, freq_all{ss}, psd_tAMM{ss}, 'Color', col, 'LineWidth', 1.2);
    plot(ax_HFC,  freq_all{ss}, psd_HFC{ss},  'Color', col, 'LineWidth', 1.2);

    legend_entries{end+1} = subj_label; %#ok<AGROW>

end

% Shared y-limits
yl = [min([ax_tAMM.YLim(1), ax_HFC.YLim(1)]), ...
      max([ax_tAMM.YLim(2), ax_HFC.YLim(2)])];
ax_tAMM.YLim = yl;
ax_HFC.YLim  = yl;

% Labels and formatting
for ax = [ax_tAMM, ax_HFC]
    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, 'Power (log_{10})');
    box(ax, 'off');
    set(ax, 'FontSize', 12);
end

title(ax_tAMM, 'tAMM');
title(ax_HFC,  'HFC');
title(t, 'Group PSD — channel average (bad channels removed)');

%legend(ax_HFC, legend_entries, 'Location', 'best', 'FontSize', 9);

%% SAVE

saveas(fig, [p.directories.base_dir '/analysis/group_psd_tAMM_vs_HFC.png']);
fprintf('Saved group PSD figure to: %s\n', p.badchannels.group_psd_comparison);
close(fig);

%%
%% SECOND FIGURE — group mean ± SD, both methods on one axis

% Stack PSDs across subjects: [freq x n_subs]
valid = ~cellfun(@isempty, psd_tAMM);
freq_plot = freq_all{find(valid, 1)};

mat_tAMM = cell2mat(psd_tAMM(valid)');   % [freq x n_subs]
mat_HFC  = cell2mat(psd_HFC(valid)');

% Mean and SD across subjects
mu_tAMM  = mean(mat_tAMM, 2);
sd_tAMM  = std(mat_tAMM, 0, 2);
mu_HFC   = mean(mat_HFC,  2);
sd_HFC   = std(mat_HFC,  0, 2);

fig2 = figure('units', 'normalized', 'outerposition', [0 0 0.6 0.55]);
fig2.Name = 'Group mean PSD: tAMM vs HFC';
ax = axes(fig2); hold on;

% Force column vectors
freq_plot = freq_plot(:);
mu_tAMM   = mu_tAMM(:);   sd_tAMM = sd_tAMM(:);
mu_HFC    = mu_HFC(:);    sd_HFC  = sd_HFC(:);

% Shaded SD regions
fill(ax, [freq_plot; flipud(freq_plot)], ...
    [mu_tAMM - sd_tAMM; flipud(mu_tAMM + sd_tAMM)], ...
    [0.2 0.5 0.9], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

fill(ax, [freq_plot; flipud(freq_plot)], ...
    [mu_HFC - sd_HFC; flipud(mu_HFC + sd_HFC)], ...
    [0.9 0.3 0.2], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Mean lines
plot(ax, freq_plot, mu_tAMM, 'Color', [0.2 0.5 0.9], 'LineWidth', 2);
plot(ax, freq_plot, mu_HFC,  'Color', [0.9 0.3 0.2], 'LineWidth', 2);

% Labels
xlabel(ax, 'Frequency (Hz)');
ylabel(ax, 'Power (log_{10})');
title(ax, sprintf('Group mean PSD ± SD  (n = %d)', sum(valid)));
legend(ax, '_', '_', 'tAMM', 'HFC', 'Location', 'best', 'FontSize', 11);
box(ax, 'off');
set(ax, 'FontSize', 12);

%%

% Save
save_dir = fullfile(p.directories.base_dir, 'analysis');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
save_path = fullfile(save_dir, 'TDC092Face_group_psd_tAMM_vs_HFC_avg.png');
saveas(fig2, save_path);
fprintf('Saved group average PSD figure to: %s\n', save_path);
close(fig2);

end


%% -----------------------------------------------------------------------
%  LOCAL HELPER: fix_grad
%  Trims all length-matched grad fields to the intersected channel indices
% -----------------------------------------------------------------------
function data = fix_grad(data, ind)
    n_orig = length(data.grad.label);
    f = fieldnames(data.grad);
    for ii = 1:length(f)
        if size(data.grad.(f{ii}), 1) == n_orig
            data.grad.(f{ii}) = data.grad.(f{ii})(ind, :);
        end
    end
end