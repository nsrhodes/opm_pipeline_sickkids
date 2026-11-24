function flip = flip_face(p, sub_num, ts_ffg_r)

% Output directory to save figures
out_dir = [p.directories.base_dir '/analysis/' p.pipeline_name '/face_flipping'];
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Indices of the right occipital lobe & fusiform
rois = [44, 46, 48, 50, 52, 54];
ffg = 56;
num_rois = length(rois);

% Get timeseries of the occipital lobes
time = get_time(p);
[ts, ~] = get_timeseries(p, rois, sub_num);

% Squeeze if needed
if ndims(ts_ffg_r) == 3
    ts_ffg = squeeze(ts_ffg_r(1, :));
else
    ts_ffg = ts_ffg_r;
end

% Find the peaks within 20-300ms
[pks, locs, ~, pr] = findpeaks(abs(ts_ffg(time > 0.02 & time < 0.24)), time(time > 0.02 & time < 0.24));

% Find the two biggest peaks
[~, ind] = sort(pr, 'descend');
locs = locs(sort(ind(1:2)));

% Find ROIs with similar peaks
flag = false(num_rois, 1);
for rr = 1:num_rois

    % Squeeze
    ts_cur = squeeze(ts(rr, :));

    % Find the peaks within 20-300ms
    [~, locs_cur, ~, pr_cur] = findpeaks(abs(ts_cur(time > 0.02 & time < 0.24)), time(time > 0.02 & time < 0.24));

    % Find the peaks within 10ms of the first FFG peak
    first = abs(locs_cur - locs(1));
    if any(first < 0.10)
        [~, first_match] = min(first);
        first_match = locs_cur(first_match);
    else
        first_match = NaN;
    end

    % Find the peak within 10ms of the second FFG peak
    second = abs(locs_cur - locs(2));
    if any(second < 0.10)
        [~, second_match] = min(second);
        second_match = locs_cur(second_match);
    else
        second_match = NaN;
    end

    % Check whether the signs are the same
    if sign(ts_ffg(time == locs(1))) == sign(ts_cur(time == first_match)) && sign(ts_ffg(time == locs(2))) == sign(ts_cur(time == second_match))
        flag(rr) = true;
    end

end
% If we find no matches, check if the non-FFG are consistent with the
% calcarine
if all(~flag)

    % Calcarine ROI
    cal = find(rois == 44);

    % Squeeze
    ts_cal = squeeze(ts(cal, :));

    % Find the peaks within 20-300ms
    [pks, locs, ~, pr] = findpeaks(abs(ts_cal(time > 0.01 & time < 0.24)), time(time > 0.01 & time < 0.24));

    % Find the two biggest peaks
    [~, ind] = sort(pr, 'descend');
    locs = locs(sort(ind(1:2)));

    % Find ROIs with similar peaks
    flag = false(length(rois), 1);
    for rr = 1:length(rois)

        % Get the current
        ts_cur = squeeze(ts(rr, :));

        % Find the peaks within 20-300ms
        [~, locs_cur, ~, pr_cur] = findpeaks(abs(ts_cur(time > 0.01 & time < 0.24)), time(time > 0.01 & time < 0.24));

        % Find the peaks within 20ms of the first calcarine peak
        first = abs(locs_cur - locs(1));
        if any(first < 0.20)
            [~, first_match] = min(first);
            first_match = locs_cur(first_match);
        else
            first_match = NaN;
        end

        % Find the peak within 10ms of the second FFG peak
        second = abs(locs_cur - locs(2));
        if any(second < 0.10)
            [~, second_match] = min(second);
            second_match = locs_cur(second_match);
        else
            second_match = NaN;
        end

        % Check whether the signs are the same
        if sign(ts_cal(time == locs(1))) == sign(ts_cur(time == first_match)) && sign(ts_cal(time == locs(2))) == sign(ts_cur(time == second_match))
            flag(rr) = true;
        end

    end

end

% Set the ROIs
rois_cur = [rois(flag), ffg];

% Get the average timeseries
ts_temp = [squeeze(ts(flag, :)); ts_ffg];
ts_avg = mean(ts_temp, 1);

% Find the peaks
[~, locs_orig, ~, pr] = findpeaks(abs(ts_avg(time > 0.01 & time < 0.24)), time(time > 0.01 & time < 0.24));
[~, ind] = sort(pr, 'descend');
locs = locs_orig(sort(ind(1:2)));
tit = 'Average';

% Check if polarity of the peaks is the same
if ~(sign(ts_ffg(time == locs(1))) ~= sign(ts_ffg(time == locs(2))))

    % Plot the fusiform instead
    ts_avg = ts_ffg;
    tit = 'Fusiform';

    % Try on just the fusiform
    [~, locs_orig, ~, pr] = findpeaks(abs(ts_ffg(time > 0.01 & time < 0.24)), time(time > 0.01 & time < 0.24));
    [~, ind] = sort(pr, 'descend');
    locs = locs_orig(sort(ind(1:2)));

    % Check if still an issue
    if ~(sign(ts_ffg(time == locs(1))) ~= sign(ts_ffg(time == locs(2))))

        % Find the next point in time
        try
            ind = find(locs_orig == locs(2));
            locs = [locs(2), locs_orig(ind+1)];
        catch
            ind = find(locs_orig == locs(1));
            locs = [locs(1), locs_orig(ind+1)];
        end

    end

end

% Plot
figure('units', 'Normalized', 'outerposition', [0, 0, 1, 0.5]);
title(p.subject(sub_num));
for rr = 1:length(rois_cur)
    subplot(1, length(rois)+2, rr);
    plot(time, ts_temp(rr, :));
    xlim([-0.1, 0.3]);
    title(rois_cur(rr));
end

% Plot the average
subplot(1, length(rois)+2, length(rois_cur)+1);
hold on;
plot(time, ts_avg);
scatter(locs(1), ts_avg(time == locs(1)), 100, 'red', 'filled');
scatter(locs(2), ts_avg(time == locs(2)), 100, 'green', 'filled');
xlim([-0.1, 0.3]);
title(tit);
hold off

% Assign flip value
if ts_ffg(time == locs(1)) < ts_ffg(time == locs(2))
    flip = 1;
    if flip
        ts_ffg = -ts_ffg;
    end
else
    flip = 0;
end

% Save figure
saveas(gcf, [out_dir '/' p.subject(sub_num) '_' p.session(sub_num) '_' p.run(sub_num) '_occipital_ts.png']);

end

